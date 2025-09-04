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
from sequence_utils import reverse_complement, extract_sequence, check_barcode, calc_softclip_lens

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
        'qname': read.query_name,
        'flag':  read.flag,
        'rname': read.reference_name,
        'pos':   read.reference_start,
        'mapq':  read.mapping_quality,
        'cigar': read.cigarstring,
        'rnext': read.next_reference_name,
        'pnext': read.next_reference_start,
        'tlen':  read.template_length,
        'seq':   read.query_sequence,
        'qual':  read.query_qualities,
        'tags':  dict(read.tags),

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
    # read reference name like exon_inclusion_var1, exon_skipping_var1
    read_ref_var = re.search(r'[^_]+$', read['rname']).group(0)
    read_ref_exon = re.sub(r'_[^_]+$', '', read['rname'])

    read_pos = read['pos']
    read_pos_end = read['pos_end']
    read_seq = read['seq']

    read_start_pass = ( exon_positions.loc[exon_positions['var_id'] == read_ref_var, 'exon_end'].iloc[0] - read_pos > args.soft_clip )
    read_end_pass = ( read_pos_end - exon_positions.loc[exon_positions['var_id'] == read_ref_var, 'length'].iloc[[0, 1]].sum() > args.soft_clip
                      if read_ref_exon == "exon_inclusion" 
                      else read_pos_end - exon_positions.loc[exon_positions['var_id'] == read_ref_var, 'length'].iloc[0] > args.soft_clip ) 

    if read_start_pass and read_end_pass:
        barcode_seq = extract_sequence(read_seq, args.barcode_up, args.barcode_down, args.max_mismatch)
        barcode_dict = { 
            'type': read_ref_exon,
            'ref_id': read_ref_var,
            'barcode': barcode_seq
        }
        if args.barcode_check:
            check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
            if not check_res:
                barcode_dict = {}
        return read, {}, barcode_dict
    else:
        return {}, read, {}

def condition1(read1_cigar: str, soft_clip: int) -> bool:
    """
    Check if the first soft clip in read1's CIGAR string is within the allowed soft clip length.
    Parameters:
        -- read1_cigar: CIGAR string of read1
        -- soft_clip: maximum allowed soft clip length
    Returns:
        -- bool: True if the condition is met, False otherwise
    """
    if re.search(r'\d+S$', read1_cigar):
        _, last_softclip = calc_softclip_lens(read1_cigar)
        return last_softclip <= soft_clip
    else:
        return True
    
def condition2(read2_cigar: str, read2_start_pos: int, read1_end_pos: int, soft_clip: int) -> bool:
    """
    Check if read2 starts after read1 ends, considering the first soft clip in read2's CIGAR string.
    Parameters:
        -- read2_cigar: CIGAR string of read2
        -- read2_start_pos: start position of read2
        -- read1_end_pos: end position of read1
        -- soft_clip: maximum allowed soft clip length
    Returns:
        -- bool: True if read2 starts after read1 ends, False otherwise
    """
    if re.match(r'^\d+S', read2_cigar):
        first_softclip, _ = calc_softclip_lens(read2_cigar)
        if first_softclip <= soft_clip:
            return read2_start_pos >= read1_end_pos
        else:
            return False
    else:
        return read2_start_pos >= read1_end_pos

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
    read2_ref = read2['rname']
    read1_pos = read1['pos']
    read2_pos = read2['pos']
    read1_pos_end = read1['pos_end']
    read2_pos_end = read2['pos_end']
    read2_seq_rc = read2['seq_rc']
    read1_cigar = read1['cigar']
    read2_cigar = read2['cigar']

    if read1_ref == read2_ref:
        # read reference name like exon_inclusion_var1, exon_skipping_var1
        read_ref_var = re.search(r'[^_]+$', read1_ref).group(0)
        read_ref_exon = re.sub(r'_[^_]+$', '', read1_ref)

        read_start_pass = ( exon_positions.loc[exon_positions['var_id'] == read_ref_var, 'exon_end'].iloc[0] - read1_pos > args.soft_clip )
        read_end_pass = ( read2_pos_end - exon_positions.loc[exon_positions['var_id'] == read_ref_var, 'length'].iloc[[0, 1]].sum() > args.soft_clip
                          if read_ref_exon == "exon_inclusion" 
                          else read2_pos_end - exon_positions.loc[exon_positions['var_id'] == read_ref_var, 'length'].iloc[0] > args.soft_clip )

        if read_start_pass and read_end_pass:
            if condition1(read1_cigar, args.soft_clip) and \
               condition2(read2_cigar, read2_pos, read1_pos_end, args.soft_clip):
                # 1. assum read1 contains variant and read2 contains barcode
                # 2. read2 is reversed complemented
                # 3. up and down sequences are in the forward strand
                barcode_seq = extract_sequence(read2_seq_rc, args.barcode_up, args.barcode_down, args.max_mismatch)
                barcode_dict = { 
                    'type': read_ref_exon,
                    'ref_id': read_ref_var,
                    'barcode': barcode_seq
                }
                if args.barcode_check:
                    check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
                    if not check_res:
                        barcode_dict = {}
                return read1, read2, {}, {}, barcode_dict

    return {}, {}, read1, read2, {}
    
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
        pysam.AlignmentFile(mapped_bam, "wb", header = bamfile.header) as mapped_fh, \
        pysam.AlignmentFile(wrong_bam, "wb", header = bamfile.header) as unmapped_fh:
    
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
                    mapped_reads, unmapped_reads, reads_barcodes = zip(*flat_results)

                    for read_dict in mapped_reads:
                        if read_dict:
                            mapped_fh.write(dict_to_segment(read_dict))

                    for read_dict in unmapped_reads:
                        if read_dict:
                            unmapped_fh.write(dict_to_segment(read_dict))

                    read_chunk = []
                    yield pl.DataFrame(reads_barcodes)

            # process any remaining reads after file ends
            if read_chunk:
                read_batches = [
                    read_chunk[i:i + batch_size] 
                    for i in range(0, len(read_chunk), batch_size)
                ]

                args_list = [ batch for batch in read_batches]
                batch_results = list(executor.map(function_processpool_se, args_list))

                flat_results = [item for batch in batch_results for item in batch]
                mapped_reads, unmapped_reads, reads_barcodes = zip(*flat_results)

                for read_dict in mapped_reads:
                    if read_dict:
                        mapped_fh.write(dict_to_segment(read_dict))

                for read_dict in unmapped_reads:
                    if read_dict:
                        unmapped_fh.write(dict_to_segment(read_dict))
                
                read_chunk = []
                yield pl.DataFrame(reads_barcodes)
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
                    mapped_reads1, mapped_reads2, unmapped_reads1, unmapped_reads2, reads_barcodes = zip(*flat_results)

                    for read1_dict, read2_dict in zip(mapped_reads1, mapped_reads2):
                        if read1_dict and read2_dict:
                            mapped_fh.write(dict_to_segment(read1_dict))
                            mapped_fh.write(dict_to_segment(read2_dict))
                    
                    for read1_dict, read2_dict in zip(unmapped_reads1, unmapped_reads2):
                        if read1_dict and read2_dict:
                            unmapped_fh.write(dict_to_segment(read1_dict))
                            unmapped_fh.write(dict_to_segment(read2_dict))
                    
                    read_chunk = []
                    yield pl.DataFrame(reads_barcodes)
            
            # process any remaining reads after file ends
            if read_chunk:
                read_batches = [
                    read_chunk[i:i + batch_size] 
                    for i in range(0, len(read_chunk), batch_size)
                ]

                args_list = [ batch for batch in read_batches]
                batch_results = list(executor.map(function_processpool_pe, args_list))

                flat_results = [item for batch in batch_results for item in batch]
                mapped_reads1, mapped_reads2, unmapped_reads1, unmapped_reads2, reads_barcodes = zip(*flat_results)

                for read1_dict, read2_dict in zip(mapped_reads1, mapped_reads2):
                    if read1_dict and read2_dict:
                        mapped_fh.write(dict_to_segment(read1_dict))
                        mapped_fh.write(dict_to_segment(read2_dict))
                
                for read1_dict, read2_dict in zip(unmapped_reads1, unmapped_reads2):
                    if read1_dict and read2_dict:
                        unmapped_fh.write(dict_to_segment(read1_dict))
                        unmapped_fh.write(dict_to_segment(read2_dict))
                
                read_chunk = []
                yield pl.DataFrame(reads_barcodes)

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process a bwa bam file for canonical splicing events.", allow_abbrev = False)
    parser.add_argument("--bam_file",            type = str, required = True,       help = "bam file sorted by read name")
    parser.add_argument("--barcode_file",        type = str, required = True,       help = "barcode association file")
    parser.add_argument("--exon_pos",            type = str, required = True,       help = "exon position file")
    parser.add_argument("--read_type",           type = str, default = 'se',        help = "sequence read type (se or pe)", choices = ['se', 'pe'],)
    parser.add_argument("--soft_clip",           type = int, default = 5,           help = "soft clip tolerance")
    parser.add_argument("--barcode_up",          type = str, required = True,       help = "sequence before barcode in read2")
    parser.add_argument("--barcode_down",        type = str, required = True,       help = "sequence after barcode in read2")
    parser.add_argument("--barcode_check",       action="store_true",               help = "enable barcode checking against template")
    parser.add_argument("--barcode_temp",        type = str,                        help = "template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNNATNNNNATNNNNATNNNN')")
    parser.add_argument("--max_mismatch",        type = int, default = 2,           help = "max mismatches allowed in up/down matches")
    parser.add_argument("--barcode_mismatch",    type = int, default = 1,           help = "number of mismatches allowed in barcode checking")
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
    exon_positions = pd.read_csv(args.exon_pos, sep = "\t", header = None)
    exon_positions.columns = ["var_id", "exon_id", "exon_start", "exon_end"]
    exon_positions["var_id"] = exon_positions["var_id"].astype(str)
    exon_positions["exon_id"] = exon_positions["exon_id"].astype(str)
    exon_positions["exon_start"] = exon_positions["exon_start"].astype(int)
    exon_positions["exon_end"] = exon_positions["exon_end"].astype(int)
    exon_positions["length"] = exon_positions["exon_end"] - exon_positions["exon_start"] + 1

    bar_var_df = pl.read_csv(args.barcode_file, separator = "\t", has_header = True, columns = ["barcode", "varid"] )

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    mapped_bam = f"{output_prefix}.filtered.bam"
    if os.path.exists(mapped_bam):
        os.remove(mapped_bam)

    wrong_bam = f"{output_prefix}.wrongmap.bam"
    if os.path.exists(wrong_bam):
        os.remove(wrong_bam)

    barcode_out = f"{output_prefix}.barcodes.txt"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)

    # -- parallel processing -- #
    barcode_list = []
    bamfile = pysam.AlignmentFile(args.bam_file, "rb")
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing bam file, please wait...", flush = True)
    for i, chunk_result in enumerate(read_bam_in_chunk(args.bam_file, args.read_type, args.chunk_size, args.threads)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {len(chunk_result)} reads", flush = True)
        barcode_list.append(chunk_result)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)

    # -- clean and format the extracted barcodes from reads -- #
    filtered_barcode_list = [df for df in barcode_list if df.shape[1] > 0]
    if filtered_barcode_list:
        barcode_df = pl.concat(filtered_barcode_list, how = "vertical")
        barcode_df = barcode_df.filter(pl.col("barcode").is_not_null())
    else:
        barcode_df = pl.DataFrame() 
    
    if barcode_df.is_empty():
        with open(barcode_out, "w") as f:
            f.write("no barcode found in the reads, please check your barcode marker or template!\n")
    else:
        barcode_varid_df = barcode_df.join(bar_var_df, on = "barcode", how = "left")
        barcode_varid_df.write_csv(barcode_out, separator = "\t")
