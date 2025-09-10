#-- import modules --#
import os
import sys
import argparse
import re
import pysam
import csv
import pandas as pd
import polars as pl
from Bio import SeqIO
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from collections import defaultdict
from sequence_utils import reverse_complement, extract_sequence, check_barcode

#-- functions --#
def extract_read_info(read: pysam.AlignedSegment) -> dict:
    """
    Parse pysam read object to extract relevant information.
    Parameters:
        -- read: a pysam AlignedSegment
    Returns:
        -- dict: a read dict
    """
    return {
        'qname':       read.query_name,
        'flag':        read.flag,
        'rname':       read.reference_name,
        'pos':         read.reference_start,
        'mapq':        read.mapping_quality,
        'cigar':       read.cigarstring,
        'cigartuples': read.cigartuples,
        'rnext':       read.next_reference_name,
        'pnext':       read.next_reference_start,
        'tlen':        read.template_length,
        'seq':         read.query_sequence,
        'qual':        read.query_qualities,
        'tags':        dict(read.tags),

        'pos_end':    read.reference_end,
        'is_reverse': read.is_reverse,
        'seq_rc':     reverse_complement(read.query_sequence),
        'qual_rc':    read.query_qualities[::-1]
    }

def process_se_read(read: dict) -> list:
    """
    Process a single-end read to extract variants and barcodes.
    Parameters:
        -- read: dict containing read information
    Returns:
        -- list: list containing read dict, unmapped read dict, and barcode dict
    """
    read_ref = read['rname']
    read_pos = read['pos']
    read_seq = read['seq']

    barcode_seq = extract_sequence(read_seq, args.barcode_up, args.barcode_down, args.max_mismatch)
    if barcode_seq in ("upstream not found", "downstream not found"):
        barcode_seq = None
        var_id = None
    else:
        if args.barcode_check:
            check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
            if check_res is None:
                barcode_seq = None
                var_id = None
            else:
                barcode_seq = check_res
                var_id = dict_bar_var.get(barcode_seq, None)
        else:
            var_id = dict_bar_var.get(barcode_seq, None)

    if var_id is not None:
        read['rname'] = var_id
    
    if args.spliced:
        read_cigar_ops, read_cigar_length = zip(*read['cigartuples'])
        if args.lib_type == "muta_exon" or args.lib_type == "muta_intron":
            getseq_id = read_ref
        else:
            getseq_id = read['rname']

        read_spliced_list = []
        for op, length in zip(read_cigar_ops, read_cigar_length):
            if op in (0, 2): # M or D
                ref_segment = ref_sequences[getseq_id][read_pos : read_pos + length]
                read_spliced_list.append(ref_segment)
                read_pos += length  
            elif op == 3: # N
                read_pos += length
            elif op in (1, 4, 5): # I, S, H
                continue

        dict_readout = { 
            'read_ref': read_ref,
            'var_id': var_id,
            'barcode': barcode_seq,
            'spliced_sequence': ''.join(read_spliced_list) if read_spliced_list else None
        }
    else:
        dict_readout = { 
            'read_ref': read_ref,
            'var_id': var_id,
            'barcode': barcode_seq
        }

    return read, dict_readout

def process_pe_read(read_pair: tuple) -> list:
    """
    Process a paired-end read to extract variants and barcodes.
    Parameters:
        -- read: a tuple of PE read dicts containing read information
    Returns:
        -- list: list containing read dict, unmapped read dict, and barcode dict
    """
    read1, read2 = read_pair

    read1_ref = read1['rname']
    read1_pos = read1['pos']
    read1_seq = read1['seq']
    read2_ref = read2['rname']
    read2_pos = read2['pos']
    read2_seq = read2['seq']

    barcode_seq = extract_sequence(read2_seq, args.barcode_up, args.barcode_down, args.max_mismatch)
    if barcode_seq in ("upstream not found", "downstream not found"):
        barcode_seq = None
        var_id = None
    else:
        if args.barcode_check:
            check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
            if check_res is None:
                barcode_seq = None
                var_id = None
            else:
                barcode_seq = check_res
                var_id = dict_bar_var.get(barcode_seq, None)
        else:
            var_id = dict_bar_var.get(barcode_seq, None)

    if var_id is not None:
        read1['rname'] = var_id
        read2['rname'] = var_id

    if args.spliced:
        read1_cigar_ops, read1_cigar_length = zip(*read1['cigartuples'])
        read2_cigar_ops, read2_cigar_length = zip(*read2['cigartuples'])
        if args.lib_type == "muta_exon" or args.lib_type == "muta_intron":
            getseq_id = read1_ref
        else:
            getseq_id = read1['rname']

        read1_spliced_list = []
        for op, length in zip(read1_cigar_ops, read1_cigar_length):
            if op in (0, 2): # M or D
                ref_segment = ref_sequences[getseq_id][read1_pos : read1_pos + length]
                read1_spliced_list.append(ref_segment)
                read1_pos += length  
            elif op == 3: # N
                read1_pos += length
            elif op in (1, 4, 5): # I, S, H
                continue
        
        read2_spliced_list = []
        for op, length in zip(read2_cigar_ops, read2_cigar_length):
            if op in (0, 2): # M or D
                ref_segment = ref_sequences[getseq_id][read2_pos : read2_pos + length]
                read2_spliced_list.append(ref_segment)
                read2_pos += length  
            elif op == 3: # N
                read2_pos += length
            elif op in (1, 4, 5): # I, S, H
                continue

        # read1_pos and read2_pos have been updated to the end of alignments
        # get the gap sequence between read1 and read2
        # normally read1 and read2 should not overlap in the middle as reads have been separated
        if read1_pos < read2['pos']: # read2_pos now is changed to the end of read2 alignment
            pe_gap_seq = ref_sequences[getseq_id][(read1_pos + 1): read2['pos']]
        else:
            pe_gap_seq = ''
        spliced_sequence = ''.join(read1_spliced_list) + pe_gap_seq + ''.join(read2_spliced_list)

        dict_readout = { 
            'read_ref': read1_ref,
            'var_id': var_id,
            'barcode': barcode_seq,
            'spliced_sequence': spliced_sequence if spliced_sequence else None
        }
    else:
        dict_readout = { 
            'read_ref': read1_ref,
            'var_id': var_id,
            'barcode': barcode_seq
        }

    return read1, read2, dict_readout
    
def batch_process_se_reads(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads: list of read dicts
    Returns:
        -- list: list of variant dicts for the batch
    """
    results = []
    for read in batch_reads:
        result = process_se_read(read)
        results.append(result)
    return results

def batch_process_pe_reads(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads: list of read dicts
    Returns:
        -- list: list of variant dicts for the batch
    """
    results = []
    for read_pair in batch_reads:
        result = process_pe_read(read_pair)
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
    return batch_process_pe_reads(args)

def dict_to_segment(read: dict) -> pysam.AlignedSegment:
    """
    Convert a read dict to a pysam AlignedSegment.
    Parameters:
        -- read: dict containing read information
    Returns:
        -- pysam.AlignedSegment object
    """
    segment = pysam.AlignedSegment(header = bamfile.header)
    segment.query_name = read['qname']
    segment.flag = read['flag']
    segment.reference_name = read['rname']
    segment.reference_start = read['pos']
    segment.mapping_quality = read['mapq']
    segment.cigarstring = read['cigar']
    segment.next_reference_name = read['rnext']
    segment.next_reference_start = read['pnext']
    segment.template_length = read['tlen']
    segment.query_sequence = read['seq']
    segment.query_qualities = read['qual']
    segment.tags = list(read['tags'].items())
    
    return segment
   
def read_bam_in_chunk(bam_file: str, read_type: str, chunk_size: int, threads: int):
    """
    Read single-end BAM file in chunks and process reads.
    Parameters:
        -- bam_file: Path to the BAM file
        -- read_type: Type of reads ('se' for single-end, 'pe' for paired-end)
        -- chunk_size: Number of reads to process in each chunk
        -- threads: Number of threads to use for processing
    Returns:
        Generator yielding processed reads.
    """
    with pysam.AlignmentFile(bam_file, "rb", threads = threads) as bam_handle, \
        ProcessPoolExecutor(max_workers=threads) as executor, \
        pysam.AlignmentFile(fixed_bam, "wb", header = bamfile.header) as fixed_fh:
    
        batch_size = min(chunk_size, 5000)
        
        # process single-end reads
        if read_type == 'se':
            read_chunk = []
            for read in bam_handle.fetch(until_eof = True):
                read_chunk.append(extract_read_info(read))

                if len(read_chunk) >= chunk_size:
                    read_batches = [
                        read_chunk[i:i + batch_size] 
                        for i in range(0, len(read_chunk), batch_size)
                    ]

                    args_list = [ batch for batch in read_batches]
                    batch_results = list(executor.map(function_processpool_se, args_list))

                    flat_results = [item for batch in batch_results for item in batch]

                    read_chunk = []
                    fixed_reads, out_barcodes = zip(*flat_results)
                    for dict_read in fixed_reads:
                        fixed_fh.write(dict_to_segment(dict_read))
                    yield pl.DataFrame(out_barcodes)

            # process any remaining reads after file ends
            if read_chunk:
                read_batches = [
                    read_chunk[i:i + batch_size] 
                    for i in range(0, len(read_chunk), batch_size)
                ]

                args_list = [ batch for batch in read_batches]
                batch_results = list(executor.map(function_processpool_se, args_list))

                flat_results = [item for batch in batch_results for item in batch]

                read_chunk = []
                fixed_reads, out_barcodes = zip(*flat_results)
                for dict_read in fixed_reads:
                    fixed_fh.write(dict_to_segment(dict_read))
                yield pl.DataFrame(out_barcodes)
        # process paired-end reads
        else:
            read_chunk = []
            current_read = None
            for read in bam_handle.fetch(until_eof = True):
                if current_read is None:
                    current_read = read
                    continue

                if read.query_name == current_read.query_name:
                    if read.is_read1:
                        read_chunk.append((extract_read_info(read), extract_read_info(current_read)))
                    else:
                        read_chunk.append((extract_read_info(current_read), extract_read_info(read)))
                    current_read = None
                else:
                    raise ValueError(f"PE read names differ, make sure bam file to be sorted by read name!")
                
                if len(read_chunk) >= chunk_size:
                    read_batches = [
                        read_chunk[i:i + batch_size] 
                        for i in range(0, len(read_chunk), batch_size)
                    ]

                    args_list = [ batch for batch in read_batches]
                    batch_results = list(executor.map(function_processpool_pe, args_list))

                    flat_results = [item for batch in batch_results for item in batch]
                    
                    read_chunk = []
                    fixed_reads_r1, fixed_reads_r2, out_barcodes = zip(*flat_results)
                    for dict_read1, dict_read2 in zip(fixed_reads_r1, fixed_reads_r2):
                        fixed_fh.write(dict_to_segment(dict_read1))
                        fixed_fh.write(dict_to_segment(dict_read2))
                    yield pl.DataFrame(out_barcodes)
            
            # process any remaining reads after file ends
            if read_chunk:
                read_batches = [
                    read_chunk[i:i + batch_size] 
                    for i in range(0, len(read_chunk), batch_size)
                ]

                args_list = [ batch for batch in read_batches]
                batch_results = list(executor.map(function_processpool_pe, args_list))

                flat_results = [item for batch in batch_results for item in batch]

                read_chunk = []
                fixed_reads_r1, fixed_reads_r2, out_barcodes = zip(*flat_results)
                for dict_read1, dict_read2 in zip(fixed_reads_r1, fixed_reads_r2):
                    fixed_fh.write(dict_to_segment(dict_read1))
                    fixed_fh.write(dict_to_segment(dict_read2))
                yield pl.DataFrame(out_barcodes)

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process a hisat2 bam file for novel splicing events.", allow_abbrev = False)
    parser.add_argument("--lib_type",            type = str, required = True,       help = "library type", choices = ['random_intron', 'random_exon', 'muta_intron', 'muta_exon'])
    parser.add_argument("--bam_file",            type = str, required = True,       help = "bam file sorted by read name")
    parser.add_argument("--ref_file",            type = str, required = True,       help = "reference fasta file")
    parser.add_argument("--barcode_file",        type = str, required = True,       help = "barcode association file")
    parser.add_argument("--read_type",           type = str, default = 'se',        help = "sequence read type (se or pe)", choices = ['se', 'pe'])
    parser.add_argument("--barcode_up",          type = str, required = True,       help = "sequence before barcode in read2")
    parser.add_argument("--barcode_down",        type = str, required = True,       help = "sequence after barcode in read2")
    parser.add_argument("--barcode_check",       action="store_true",               help = "enable barcode checking against template")
    parser.add_argument("--barcode_temp",        type = str,                        help = "template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNNATNNNNATNNNNATNNNN')")
    parser.add_argument("--max_mismatch",        type = int, default = 2,           help = "max mismatches allowed in up/down matches")
    parser.add_argument("--barcode_mismatch",    type = int, default = 1,           help = "number of mismatches allowed in barcode checking")
    parser.add_argument("--spliced",             action="store_true",               help = "enable to create spliced products/seqeuences")
    parser.add_argument("--output_dir",          type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",       type = str, default = '',          help = "output prefix")
    parser.add_argument("--chunk_size",          type = int, default = 1000,        help = "chunk size for processing reads")
    parser.add_argument("--threads",             type = int, default = 4,           help = "number of threads")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    
    if args.barcode_check and args.barcode_temp is None:
        parser.error("--barcode_temp is required when --barcode_check is set")

    if args.output_prefix == '':
        output_prefix = os.path.splitext(os.path.basename(args.bam_file))[0]
    else:
        output_prefix = args.output_prefix

    # -- read input files -- #
    ref_sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.ref_file, "fasta")}

    df_bar_var = pl.read_csv(args.barcode_file, separator = "\t", has_header = True, columns = ["barcode", "var_id"])
    dict_bar_var = dict(zip(df_bar_var["barcode"], df_bar_var["var_id"]))

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    fixed_bam = f"{output_prefix}.fixed.bam"
    if os.path.exists(fixed_bam):
        os.remove(fixed_bam)

    barcode_out = f"{output_prefix}.barcodes.txt"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)

    # -- parallel processing -- #
    list_result = []
    bamfile = pysam.AlignmentFile(args.bam_file, "rb")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing bam file, please wait...", flush = True)
    for i, chunk_result in enumerate(read_bam_in_chunk(args.bam_file, args.read_type, args.chunk_size, args.threads)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {len(chunk_result)} reads", flush = True)
        list_result.append(chunk_result)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)    

    list_result_filtered = [df for df in list_result if df.shape[1] > 0]
    if list_result_filtered:
        df_result = pl.concat(list_result_filtered, how = "vertical")
        if args.spliced:
            df_result_count = df_result.group_by(["read_ref", "var_id", "barcode", "spliced_sequence"]).agg(pl.count().alias("count"))
            df_result_count.write_csv(barcode_out, separator = "\t", null_value = "NA")
        else:
            df_result_count = df_result.group_by(["read_ref", "var_id", "barcode"]).agg(pl.count().alias("count"))
            df_result_count.write_csv(barcode_out, separator = "\t", null_value = "NA")
    else:
        with open(barcode_out, "w") as f:
            f.write("no results, please check the input files and script options!\n")
