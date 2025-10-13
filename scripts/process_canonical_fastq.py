#-- import modules --#
import io
import os
import sys
import argparse
import re
import gc
import subprocess
import polars as pl
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice
from collections import Counter

from sequence_utils import (
    pigz_open,
    fastq_iter_se,
    fastq_iter_pe,
    read_first_fasta_seq,
    match_approximate,
    extract_sequence,
    reverse_complement,
    check_barcode
)

def init_worker():
    """
    Initilize worker
    """
    global global_dict_bar_var, global_pre_exon, global_post_exon
    global_dict_bar_var = dict_bar_var
    global_pre_exon = pre_exon
    global_post_exon = post_exon

def process_se_read(read: tuple) -> list:
    """
    Detect canonical splicing event in a read
    Parameters:
        -- read: tuple (header, sequence, quality) for read
    Returns:
        -- tuple: (read, barcode) 
    """
    read_seq = read[1]

    barcode_seq = extract_sequence(read_seq, args.barcode_up, args.barcode_down, args.max_mismatch)
    if (barcode_seq in {"upstream not found", "downstream not found"}):
        read_seq_rc = reverse_complement(read_seq)
        barcode_seq = extract_sequence(read_seq_rc, args.barcode_up, args.barcode_down, args.max_mismatch)
        if (barcode_seq in {"upstream not found", "downstream not found"}):
            return read, {}
        else:
            read_seq = read_seq_rc

    if args.barcode_check:
        check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
        if check_res is None:
            return read, {}
        else:
            barcode_seq = check_res

    res_bar_var = global_dict_bar_var.get(barcode_seq)    
    if res_bar_var is not None:
        exon_seq = res_bar_var["exon"]
        ref_canonical_inclusion = global_pre_exon + exon_seq + global_post_exon
        ref_canonical_skipping = global_pre_exon + global_post_exon
        ref_found = match_approximate(read_seq, ref_canonical_inclusion, args.max_mismatch, "hamming")
        if ref_found == -1:
            ref_found = match_approximate(read_seq, ref_canonical_skipping, args.max_mismatch, "hamming")
            if ref_found == -1:
                return read, {}
            else:
                dict_barcode = { 'read_ref' : res_bar_var["var_id"],
                                 'var_id'   : res_bar_var["var_id"],
                                 'barcode'  : barcode_seq,
                                 'ref_type' : "canonical_skipping" }
                return (), dict_barcode
        else:
            dict_barcode = { 'read_ref' : res_bar_var["var_id"],
                             'var_id'   : res_bar_var["var_id"],
                             'barcode'  : barcode_seq,
                             'ref_type' : "canonical_inclusion" }
            return (), dict_barcode
    else:
        return read, {}
    
def process_pe_pair(read_pair: tuple) -> list:
    """
    Detect canonical splicing event in a pair of reads
    Parameters:
        -- read: tuple (read1, read2) and read1, read2 are tuple (header, sequence, quality)
    Returns:
        -- tuple: tuple: (read1, read2, barcode) 
    """
    read1, read2 = read_pair

    read1_seq = read1[1]
    read2_seq = read2[1]

    barcode_seq = extract_sequence(read2_seq, args.barcode_up, args.barcode_down, args.max_mismatch)
    if (barcode_seq in {"upstream not found", "downstream not found"}):
        read2_seq_rc = reverse_complement(read2_seq)
        barcode_seq = extract_sequence(read2_seq_rc, args.barcode_up, args.barcode_down, args.max_mismatch)
        if (barcode_seq in {"upstream not found", "downstream not found"}):
            return read1, read2, {}
        else:
            read2_seq = read2_seq_rc
    
    if args.barcode_check:
        check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
        if check_res is None:
            return read1, read2, {}
        else:
            barcode_seq = check_res

    res_bar_var = global_dict_bar_var.get(barcode_seq)    
    if res_bar_var is not None:
        exon_seq = res_bar_var["exon"]
        ref_canonical_inclusion = global_pre_exon + exon_seq + global_post_exon
        ref_canonical_skipping = global_pre_exon + global_post_exon
        # Note: assuming read1_seq + read2_seq is the whole sequence
        # hamming distance requires the same length
        # need levenshtein distance
        ref_found = match_approximate(read1_seq + read2_seq, ref_canonical_inclusion, 4 * args.max_mismatch, "levenshtein")
        if ref_found == -1:
            ref_found = match_approximate(read1_seq + read2_seq, ref_canonical_skipping, 4 * args.max_mismatch, "levenshtein")
            if ref_found == -1:
                return read1, read2, {}
            else:
                dict_barcode = { 'read_ref' : res_bar_var["var_id"],
                                 'var_id'   : res_bar_var["var_id"],
                                 'barcode'  : barcode_seq,
                                 'ref_type' : "canonical_skipping" }
                return (), (), dict_barcode
        else:
            dict_barcode = { 'read_ref' : res_bar_var["var_id"],
                             'var_id'   : res_bar_var["var_id"],
                             'barcode'  : barcode_seq,
                             'ref_type' : "canonical_inclusion" }
            return (), (), dict_barcode
    else:
        return read1, read2, {}

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

def process_se_reads_in_chunk(path_read, fh_fastq):
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

            barcode_counts = Counter()
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
                for read, dict_barcode in batch_result:
                    if read:
                        fh_fastq.write(f"{read[0]}\n{read[1]}\n+\n{read[2]}\n".encode("utf-8"))
                    elif dict_barcode:
                        key = (dict_barcode["read_ref"], 
                               dict_barcode["var_id"],
                               dict_barcode["barcode"], 
                               dict_barcode["ref_type"])
                        barcode_counts[key] += 1
                # -- free memory -- #
                del batch_result
                gc.collect()

            if barcode_counts:
                df_yield = pl.DataFrame({"read_ref": [k[0] for k in barcode_counts.keys()],
                                         "var_id":   [k[1] for k in barcode_counts.keys()],
                                         "barcode":  [k[2] for k in barcode_counts.keys()],
                                         "ref_type": [k[3] for k in barcode_counts.keys()],
                                         "count":    [v for v in barcode_counts.values()]})
            else:
                df_yield = pl.DataFrame(schema={"read_ref": pl.Utf8, "var_id": pl.Utf8, "barcode": pl.Utf8, "ref_type": pl.Utf8, "count": pl.Int64})

            # -- free memory -- #
            del read_chunk, read_batches, futures, barcode_counts
            gc.collect()

            yield df_yield
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

            barcode_counts = Counter()
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
                for read1, read2, dict_barcode in batch_result:
                    if read1 and read2:
                        fh_fastq_r1.write(f"{read1[0]}\n{read1[1]}\n+\n{read1[2]}\n".encode("utf-8"))
                        fh_fastq_r2.write(f"{read2[0]}\n{read2[1]}\n+\n{read2[2]}\n".encode("utf-8"))
                    elif dict_barcode:
                        key = (dict_barcode["read_ref"], 
                               dict_barcode["var_id"],
                               dict_barcode["barcode"], 
                               dict_barcode["ref_type"])
                        barcode_counts[key] += 1
                # -- free memory -- #
                del batch_result
                gc.collect()

            if barcode_counts:
                df_yield = pl.DataFrame({"read_ref": [k[0] for k in barcode_counts.keys()],
                                         "var_id":   [k[1] for k in barcode_counts.keys()],
                                         "barcode":  [k[2] for k in barcode_counts.keys()],
                                         "ref_type": [k[3] for k in barcode_counts.keys()],
                                         "count":    [v for v in barcode_counts.values()]})
            else:
                df_yield = pl.DataFrame(schema={"read_ref": pl.Utf8, "var_id": pl.Utf8, "barcode": pl.Utf8, "ref_type": pl.Utf8, "count": pl.Int64})

            # -- free memory -- #
            del read_chunk, read_batches, futures, barcode_counts
            gc.collect()

            yield df_yield
    fh_read1.close()
    fh_read2.close()

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process a bwa bam file for canonical splicing events.", allow_abbrev = False)
    parser.add_argument("--reads",               type = str, required = True,       help = "fastq file(s), eg: se_reads.fq.gz or pe_r1.fq.gz,pe_r2.fq.gz")
    parser.add_argument("--read_type",           type = str, default = 'se',        help = "sequence read type (se or pe)", choices = ['se', 'pe'])
    parser.add_argument("--ref_file",            type = str, required = True,       help = "reference fasta file (reads must cover the whole reference sequence, end to end)")
    parser.add_argument("--barcode_file",        type = str, required = True,       help = "barcode association file")
    parser.add_argument("--barcode_up",          type = str, required = True,       help = "sequence before barcode in read2")
    parser.add_argument("--barcode_down",        type = str, required = True,       help = "sequence after barcode in read2")
    parser.add_argument("--barcode_check",       action="store_true",               help = "enable barcode checking against template")
    parser.add_argument("--barcode_temp",        type = str,                        help = "template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNNATNNNNATNNNNATNNNN')")
    parser.add_argument("--max_mismatch",        type = int, default = 2,           help = "max mismatches allowed in up/down matches")
    parser.add_argument("--barcode_mismatch",    type = int, default = 1,           help = "number of mismatches allowed in barcode checking")
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
    
    if args.barcode_check and args.barcode_temp is None:
        parser.error("--barcode_temp is required when --barcode_check is set")

    if args.output_prefix == '':
        output_prefix = os.path.splitext(os.path.basename(input_path_parts[0]))[0]
    else:
        output_prefix = args.output_prefix

    # -- read input files -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Reading files, please wait...", flush = True)
    df_bar_var = pl.read_csv(args.barcode_file, separator = "\t", has_header = True, columns = ["barcode", "variant", "var_id"])
    df_bar_var = df_bar_var.with_columns(pl.col("variant")
                                           .map_elements(lambda s: "".join(c for c in s if c.isupper()), return_dtype = pl.Utf8)
                                           .alias("exon"))
    df_bar_var = df_bar_var.drop("variant")
    dict_bar_var = { row["barcode"]: {"var_id": row["var_id"], "exon": row["exon"]}
                     for row in df_bar_var.iter_rows(named = True) }

    # -- free memory -- #
    del df_bar_var
    gc.collect()

    first_fasta_seq = read_first_fasta_seq(args.ref_file)
    exons_in_first_fasta_seq = re.findall(r"[A-Z]+", first_fasta_seq)
    exons_in_first_fasta_seq = [e.split("N", 1)[0] for e in exons_in_first_fasta_seq]
    pre_exon = exons_in_first_fasta_seq[0]
    post_exon = exons_in_first_fasta_seq[2]

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    if args.read_type == 'se':
        fail_fastq = f"{output_prefix}.filter_se.fail.fastq.gz"
        if os.path.exists(fail_fastq):
            os.remove(fail_fastq)
    else:
        fail_fastq_r1 = f"{output_prefix}.filter_se.fail.r1.fastq.gz"
        if os.path.exists(fail_fastq_r1):
            os.remove(fail_fastq_r1)
        fail_fastq_r2 = f"{output_prefix}.filter_se.fail.r2.fastq.gz"
        if os.path.exists(fail_fastq_r2):
            os.remove(fail_fastq_r2)

    barcode_out = f"{output_prefix}.canonical_barcodes.tsv"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)
    
    # -- parallel processing -- #
    list_results = []
    if args.read_type == 'se':
        pigz_proc = subprocess.Popen( ["pigz", "-c", "-p", str(args.threads)],
                                      stdin=subprocess.PIPE,
                                      stdout=open(fail_fastq, "wb") )
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing fastq file, please wait...", flush = True)
        with pigz_proc.stdin as fh_fastq:
            for i, chunk_result in enumerate(process_se_reads_in_chunk(path_read, fh_fastq)):
                print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {args.chunk_size} reads", flush = True)                
                if not chunk_result.is_empty():
                    list_results.append(chunk_result)
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)
        pigz_proc.wait()
    else:
        pigz_proc1 = subprocess.Popen( ["pigz", "-c", "-p", str(args.threads)],
                                       stdin=subprocess.PIPE,
                                       stdout=open(fail_fastq_r1, "wb") )
        pigz_proc2 = subprocess.Popen( ["pigz", "-c", "-p", str(args.threads)],
                                       stdin=subprocess.PIPE,
                                       stdout=open(fail_fastq_r2, "wb") )
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing fastq file, please wait...", flush = True)
        with pigz_proc1.stdin as fh_fastq_r1, pigz_proc2.stdin as fh_fastq_r2:
            for i, chunk_result in enumerate(process_pe_pairs_in_chunk(path_read1, path_read2, fh_fastq_r1, fh_fastq_r2)):
                print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {args.chunk_size} reads", flush = True)                
                if not chunk_result.is_empty():
                    list_results.append(chunk_result)
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)
        pigz_proc1.wait()
        pigz_proc2.wait()

    # -- clean and format the extracted barcodes from reads -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating barcode results, please wait...", flush = True)
    list_results_filtered = [df for df in list_results if df.height > 0]
    if list_results_filtered:
        df_barcode = pl.concat(list_results_filtered, how = "vertical")
        df_barcode_counts = ( df_barcode.group_by(["read_ref", "var_id", "barcode", "ref_type"])
                                        .agg(pl.sum("count").alias("count"))
                                        .select(["read_ref", "var_id", "barcode", "count", "ref_type"]) )
        df_barcode_counts.write_csv(barcode_out, separator = "\t", null_value = "NA")
    else:
        with open(barcode_out, "w") as f:
            f.write("no barcode found in the reads, please check your barcode marker or template!\n")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done.", flush = True)
