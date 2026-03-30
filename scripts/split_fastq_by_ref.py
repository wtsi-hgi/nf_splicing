#-- import modules --#
import gzip
import io
import os
import sys
import argparse
import re
import gc
import subprocess
from Bio import SeqIO
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice

from sequence_utils import (
    pigz_open,
    fastq_iter_se,
    fastq_iter_pe
)

def init_worker():
    """
    Initilize worker
    """
    global global_dict_ref
    global_dict_ref = dict_ref

def process_se_read(read: tuple) -> list:
    """
    Detect canonical splicing event in a read
    Parameters:
        -- read: tuple (header, sequence, quality) for read
    Returns:
        -- tuple: (read, barcode) 
    """
    var_id = read[0].split("||")[-1].strip()
    exon_id = "_".join(var_id.split("_")[:2])

    if "_del" in var_id:
        m = re.search(r"_del(\d+)to(\d+)", var_id)
        if m:
            start = int(m.group(1))
            end = int(m.group(2))
            if (end - start + 1) == 21:
                return var_id, read

    return exon_id, read
    
def process_pe_pair(read_pair: tuple) -> list:
    """
    Detect canonical splicing event in a pair of reads
    Parameters:
        -- read: tuple (read1, read2) and read1, read2 are tuple (header, sequence, quality)
    Returns:
        -- tuple: tuple: (read1, read2, barcode) 
    """
    read1, read2 = read_pair
    var_id = read1[0].split("||")[-1].strip()
    exon_id = "_".join(var_id.split("_")[:2])

    if "_del" in var_id:
        m = re.search(r"_del(\d+)to(\d+)", var_id)
        if m:
            start = int(m.group(1))
            end = int(m.group(2))
            if (end - start + 1) == 21:
                return var_id, read1, read2

    return exon_id, read1, read2

def batch_process_se_reads(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads
    Returns:
        -- list of tuples
    """
    results = []
    for read in batch_reads:
        result = process_se_read(read)
        results.append(result)
    return results

def batch_process_pe_pairs(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads
    Returns:
        -- list of tuples
    """
    results = []
    for read_pair in batch_reads:
        result = process_pe_pair(read_pair)
        results.append(result)
    return results

def function_processpool_se(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_se_reads(args)

def function_processpool_pe(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_pe_pairs(args)

def process_se_reads_in_chunk(path_read, handles):
    """
    Read single-end FASTQ files in chunks and process reads in parallel
    Parameters:
        -- path_read: path to single-end read FASTQ file
    Yields:
        -- DataFrame: barcode
    """
    fh_read = io.TextIOWrapper(pigz_open(path_read).stdout) if path_read.endswith(".gz") else open(path_read)
    read_iter = fastq_iter_se(fh_read)

    with ProcessPoolExecutor(max_workers = args.threads, initializer = init_worker) as executor:
        while True:
            read_chunk = list(islice(read_iter, args.chunk_size))
            if not read_chunk:
                break

            # Divide chunk into batches
            # if process_long_read is very fast, we can use a larger batch size to make better use of CPU resources
            # if process_long_read is very slow, we can use a smaller batch size to make better use of CPU resources
            batch_size = min(args.chunk_size, 5000)
            read_batches = [
                read_chunk[i:i+batch_size]
                for i in range(0, len(read_chunk), batch_size)
            ]

            # executor.map() --> memeory may grow accumulatively
            # 1. creates all tasks upfront for the entire input list (read_batches)
            # 2. stores them all in an internal queue inside the Executor.
            # 3. waits for results in order — not as they finish.
            # 
            # as_completed()
            # 1. Yields each Future as soon as it finishes, regardless of submission order
            # 2. You can process and discard the result immediately.
            futures = [ executor.submit(function_processpool_se, batch) for batch in read_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                for rid, read in batch_result:
                    if rid in handles:
                        handle = handles[rid]
                        handle.write(f"{read[0]}\n{read[1]}\n+\n{read[2]}\n")

                # -- free memory -- #
                del batch_result
                gc.collect()

            yield []
    fh_read.close()

def process_pe_pairs_in_chunk(path_read1, path_read2, fh_fastq_r1, fh_fastq_r2):
    """
    Read paired-end FASTQ files in chunks and process reads in parallel
    Parameters:
        -- path_read1: path to paired-end read1 FASTQ file
        -- path_read2: path to paired-end read2 FASTQ file
    Yields:
        -- DataFrame: barcode
    """
    fh_read1 = io.TextIOWrapper(pigz_open(path_read1).stdout) if path_read1.endswith(".gz") else open(path_read1)
    fh_read2 = io.TextIOWrapper(pigz_open(path_read2).stdout) if path_read2.endswith(".gz") else open(path_read2)
    read_iter = fastq_iter_pe(fh_read1, fh_read2)

    with ProcessPoolExecutor(max_workers = args.threads, initializer = init_worker) as executor:
        while True:
            read_chunk = list(islice(read_iter, args.chunk_size))
            if not read_chunk:
                break

            # Divide chunk into batches
            # if process_long_read is very fast, we can use a larger batch size to make better use of CPU resources
            # if process_long_read is very slow, we can use a smaller batch size to make better use of CPU resources
            batch_size = min(args.chunk_size, 5000)
            read_batches = [
                read_chunk[i:i+batch_size]
                for i in range(0, len(read_chunk), batch_size)
            ]

            # executor.map() --> memeory may grow accumulatively
            # 1. creates all tasks upfront for the entire input list (read_batches)
            # 2. stores them all in an internal queue inside the Executor.
            # 3. waits for results in order — not as they finish.
            # 
            # as_completed()
            # 1. Yields each Future as soon as it finishes, regardless of submission order
            # 2. You can process and discard the result immediately.
            futures = [ executor.submit(function_processpool_pe, batch) for batch in read_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                for rid, read1, read2 in batch_result:
                    if rid in fh_fastq_r1 and rid in fh_fastq_r2:
                        handle_r1 = fh_fastq_r1[rid]
                        handle_r1.write(f"{read1[0]}\n{read1[1]}\n+\n{read1[2]}\n")
                        handle_r2 = fh_fastq_r2[rid]
                        handle_r2.write(f"{read2[0]}\n{read2[1]}\n+\n{read2[2]}\n")

                # -- free memory -- #
                del batch_result
                gc.collect()

            yield []
    fh_read1.close()
    fh_read2.close()

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process a bwa bam file for canonical splicing events.", allow_abbrev = False)
    parser.add_argument("--reads",               type = str, required = True,       help = "fastq file(s), eg: se_reads.fq.gz or pe_r1.fq.gz,pe_r2.fq.gz")
    parser.add_argument("--read_type",           type = str, default = 'se',        help = "sequence read type (se or pe)", choices = ['se', 'pe'])
    parser.add_argument("--ref_file",            type = str, required = True,       help = "reference fasta file (reads must cover the whole reference sequence, end to end)")
    parser.add_argument("--output_dir",          type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",       type = str, default = '',          help = "output prefix")
    parser.add_argument("--chunk_size",          type = int, default = 100000,      help = "chunk size for processing reads")
    parser.add_argument("--threads",             type = int, default = 40,          help = "number of threads")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    input_path = args.reads.strip()
    input_path_parts = [r.strip() for r in input_path.split(",") if r.strip()]
    if args.read_type == "se":
        if len(input_path_parts) != 1:
            raise ValueError(f"Single-end mode (se) expects 1 read file, but got {len(input_path_parts)}: {input_path}")
        path_read = input_path_parts[0]
    elif args.read_type == "pe":
        if len(input_path_parts) != 2:
            raise ValueError(f"Paired-end mode (pe) expects 2 comma-separated files, but got {len(input_path_parts)}: {input_path}")
        path_read1, path_read2 = input_path_parts
    else:
        raise ValueError(f"Invalid read type: {args.read_type} (expected 'se' or 'pe')")

    if args.output_prefix == '':
        output_prefix = os.path.splitext(os.path.basename(input_path_parts[0]))[0]
    else:
        output_prefix = args.output_prefix

    # -- read input files -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Reading files, please wait...", flush = True)
    dict_ref = {}
    for record in SeqIO.parse(args.ref_file, "fasta"):
        if "_wt" in record.id:
            dict_ref[record.id.split("_wt")[0]] = str(record.seq)
        elif "_del" in record.id:
            m = re.search(r"_del(\d+)to(\d+)", record.id)
            if m:
                start = int(m.group(1))
                end = int(m.group(2))
                if (end - start + 1) == 21:
                    dict_ref[record.id] = str(record.seq)

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    for id, seq in dict_ref.items():
        out_file = f"{output_prefix}.{id}.exon.fasta"
        with open(out_file, "w") as fh:
            fh.write(f">{id}_wt\n{seq}\n")

    if args.read_type == 'se':
        handles = {}
        for id in dict_ref.keys():
            out = f"{output_prefix}.{id}.exon.fastq.gz"
            if os.path.exists(out):
                os.remove(out)
            handles[id] = gzip.open(out, "wt")
    else:
        handles_r1 = {}
        handles_r2 = {}
        for id in dict_ref.keys():
            out_r1 = f"{output_prefix}.{id}.r1.exon.fastq.gz"
            out_r2 = f"{output_prefix}.{id}.r2.exon.fastq.gz"
            if os.path.exists(out_r1):
                os.remove(out_r1)
            if os.path.exists(out_r2):
                os.remove(out_r2)
            handles_r1[id] = gzip.open(out_r1, "wt")
            handles_r2[id] = gzip.open(out_r2, "wt")

    # -- parallel processing -- #
    if args.read_type == 'se':
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing fastq file, please wait...", flush = True)
        for i, _ in enumerate(process_se_reads_in_chunk(path_read, handles)):
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {args.chunk_size} reads", flush = True)                
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)
    else:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing fastq file, please wait...", flush = True)
        for i, _ in enumerate(process_pe_pairs_in_chunk(path_read1, path_read2, handles_r1, handles_r2)):
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {args.chunk_size} reads", flush = True)                
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done.", flush = True)
