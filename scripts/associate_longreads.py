#-- import modules --#
import io
import os
import sys
import argparse
import re
import gc
import subprocess
import numpy as np
import polars as pl
import Levenshtein
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice

#-------------------------------------------------------
# sequence processing functions
#-------------------------------------------------------
def pigz_open(path: str):
    """
    Open a gzip file using pigz for faster decompression
    Parameters:
        -- path: file path to the gzip file
    Returns:
        -- io_wrapper: a TextIOWrapper for reading the decompressed file
    """
    return subprocess.Popen(["pigz", "-dc", path], stdout = subprocess.PIPE)

def fastq_iter_se(handle):
    """
    FASTQ single-end parser yielding (header, seq, qual)
    Parameters:
        -- handle: file handle for the FASTQ file
    Yields:
        -- (header, seq, qual): a tuple containing the header, sequence, and quality
    """
    while True:
        header = handle.readline()
        if not header:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not (seq and plus and qual):
            raise ValueError("Truncated FASTQ record")

        yield (header.rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n"))

def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.
    Parameters:
        -- seq: the DNA sequence to reverse complement
    Returns:
        -- str: the reverse complement of the sequence
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def match_hamming_numpy(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq by hamming distance allowing max_mismatches
    Parameters:
        -- seq: the target sequence
        -- pattern: the pattern to match
    Returns:
        -- int: the Hamming distance, or the maximum length if they differ in length
    """
    k = len(pattern)
    n = len(seq)
    if k > n:
        return -1

    # convert the sequence to np.uint8 arrays of ASCII codes
    seq_arr = np.frombuffer(seq.encode("ascii"), dtype = np.uint8)
    pat_arr = np.frombuffer(pattern.encode("ascii"), dtype = np.uint8)

    # sliding window: create a 2D view of seq_arr of shape (n-k+1, k)
    windows = np.lib.stride_tricks.sliding_window_view(seq_arr, window_shape = k)

    # calculate hamming distances vectorized
    mismatches = np.sum(windows != pat_arr, axis = 1)
    matches = np.where(mismatches <= max_mismatches)[0]
    return int(matches[0]) if matches.size > 0 else -1

def match_levenshtein(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq by levenshtein distance allowing max_mismatches
    Parameters:
        -- seq: the sequence to search in
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: start index of the match or -1 if not found
    """
    k = len(pattern)
    n = len(seq)
    if k > n:
        return -1

    for i in range(n - k + 1):
        window = seq[i:i+k]
        if Levenshtein.distance(window, pattern) <= max_mismatches:
            return i
    return -1

def match_approximate(seq: str, pattern: str, max_mismatches: int, distance: str) -> int:
    """
    Hybrid approximate match supporting Hamming and Levenshtein
    Parameters:
        -- seq: the sequence to search in
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: start index of the match or -1 if not found
    """
    if distance == "hamming":
        return match_hamming_numpy(seq, pattern, max_mismatches)
    elif distance == "levenshtein":
        return match_levenshtein(seq, pattern, max_mismatches)
    else:
        raise ValueError(f"Unknown distance metric: {distance}")

def extract_sequence(seq: str, up_seq: str, down_seq: str, max_mismatches: int) -> str:
    """
    Extract substring from seq between approximate matches of up_seq and down_seq.
    Parameters:
        -- seq: the sequence to search in
        -- up_seq: upstream sequence to match
        -- down_seq: downstream sequence to match
        -- max_mismatches: maximum number of mismatches allowed for both matches
    Returns:
        -- str: the extracted sequence or an error message if not found
    """
    start_idx = match_approximate(seq, up_seq, max_mismatches, "hamming")
    if start_idx == -1:
        return "upstream not found"
    start_idx += len(up_seq)

    end_idx = match_approximate(seq[start_idx:], down_seq, max_mismatches, "hamming")
    if end_idx == -1:
        return "downstream not found"
    end_idx += start_idx

    return seq[start_idx:end_idx]

def check_barcode(barcode_seq: str, barcode_temp: str, max_mismatches: int):
    """
    Check if barcode_seq matches the barcode_temp allowing max_mismatches
    Parameters:
        -- barcode_seq: the sequence to check
        -- barcode_temp: the template sequence to match against
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- str or None: corrected barcode sequence if matches, else None
    """
    if len(barcode_seq) != len(barcode_temp):
        return None

    mismatch_count = 0
    barcode_corrected = []

    for s_char, p_char in zip(barcode_seq, barcode_temp):
        if p_char == 'N':
            barcode_corrected.append(s_char)
        elif s_char != p_char:
            mismatch_count += 1
            barcode_corrected.append(p_char)
            if mismatch_count > max_mismatches:
                return None
        else:
            barcode_corrected.append(s_char)

    return "".join(barcode_corrected)

#-------------------------------------------------------
# parallel processing functions
#-------------------------------------------------------
def process_long_read(long_read: tuple) -> list:
    """
    extract sequence by up/down markers in a long read
    Parameters:
        -- read: tuple (header, sequence, quality) for read
    Returns:
        -- tuple: (barcode_seq, vhh_seq, target_seq) 
    """
    long_read = long_read[1]

    barcode_seq = extract_sequence(long_read, args.barcode_up, args.barcode_down, args.max_mismatches)
    vhh_seq     = extract_sequence(long_read, args.vhh_up,     args.vhh_down,     args.max_mismatches)
    target_seq  = extract_sequence(long_read, args.target_up,  args.target_down,  args.max_mismatches)

    not_found = ["upstream not found", "downstream not found"]
    flags = [ barcode_seq in not_found, vhh_seq in not_found, target_seq in not_found ]
    if sum(flags) >= 2:
        long_read_rc = reverse_complement(long_read)
        barcode_seq = extract_sequence(long_read_rc, args.barcode_up, args.barcode_down, args.max_mismatches)
        vhh_seq     = extract_sequence(long_read_rc, args.vhh_up,     args.vhh_down,     args.max_mismatches)
        target_seq  = extract_sequence(long_read_rc, args.target_up,  args.target_down,  args.max_mismatches)
        
    return (barcode_seq, vhh_seq, target_seq)

def batch_process_long_reads(batch_long_reads: list) -> list:
    """
    Process a batch of long reads to extract sequences.
    Parameters:
        -- batch_long_reads
    Returns:
        -- list of tuples
    """
    results = []
    for long_read in batch_long_reads:
        result = process_long_read(long_read)
        results.append(result)
    return results

def function_processpool_long(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_long_reads(args)

def process_long_reads_in_chunk(path_long_reads: str):
    """
    Read long read FASTQ file in chunks and process reads in parallel
    Parameters:
        -- path_long_reads: path to long read FASTQ file
    Yields:
        -- DataFrame: barcode
    """
    fh_long_reads = io.TextIOWrapper(pigz_open(path_long_reads).stdout) if path_long_reads.endswith(".gz") else open(path_long_reads)
    read_iter = fastq_iter_se(fh_long_reads)

    with ProcessPoolExecutor(max_workers = args.threads) as executor:
        while True:
            read_chunk = list(islice(read_iter, args.chunk_size))
            if not read_chunk:
                break

            batch_size = min(args.chunk_size, 5000)
            read_batches = [
                read_chunk[i:i+batch_size]
                for i in range(0, len(read_chunk), batch_size)
            ]

            list_barcodes = []
            futures = [ executor.submit(function_processpool_long, batch) for batch in read_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                list_barcodes.append(pl.DataFrame(batch_result, schema = ["barcode_seq", "vhh_seq", "target_seq"], orient = "row"))

                # -- free memory -- #
                del batch_result
                gc.collect()

            df_yield = ( pl.concat(list_barcodes, how = "vertical")
                           .group_by(["barcode_seq", "vhh_seq", "target_seq"])
                           .agg(pl.len().alias("count")) )

            # -- free memory -- #
            del read_chunk, read_batches, futures, list_barcodes
            gc.collect()

            yield df_yield
    fh_long_reads.close()

#-------------------------------------------------------
# main function
#-------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Association using long reads.", allow_abbrev = False)
    parser.add_argument("--longread",         type = str,       required = True,       help = "the path of long read FASTQ file")
    parser.add_argument("--barcode_up",       type = str.upper, required = True,       help = "Upstream flank sequence of barcode in long read")
    parser.add_argument("--barcode_down",     type = str.upper, required = True,       help = "Downstream flank sequence of barcode in long read")
    parser.add_argument("--barcode_temp",     type = str.upper, default = None,        help = "Template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNN')")
    parser.add_argument("--barcode_mismatch", type = int,       default = 1,           help = "Number of mismatches allowed in barcode checking")
    parser.add_argument("--vhh_up",           type = str.upper, required = True,       help = "Upstream flank sequence of vhh in long read")
    parser.add_argument("--vhh_down",         type = str.upper, required = True,       help = "Downstream flank sequence of vhh in long read")
    parser.add_argument("--vhh_len",          type = str,       default = None,        help = "Length of the VNN sequence, eg., '16' or '16,17,18'")
    parser.add_argument("--target_up",        type = str.upper, required = True,       help = "Upstream flank sequence of target in long read")
    parser.add_argument("--target_down",      type = str.upper, required = True,       help = "Downstream flank sequence of target in long read")
    parser.add_argument("--target_len",       type = str,       default = None,        help = "Length of the target sequence, eg., '100' or '100,101,102'")    
    parser.add_argument("--max_mismatches",   type = int,       default = 2,           help = "Max mismatches allowed in up/down matches")
    parser.add_argument("--min_barcov",       type = int,       default = 2,           help = "Minimum coverage for association")
    parser.add_argument("--output_dir",       type = str,       default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",    type = str,       required = True,       help = "output prefix")
    parser.add_argument("--chunk_size",       type = int,       default = 100000,      help = "Chunk size for processing reads")
    parser.add_argument("--threads",          type = int,       default = 40,          help = "Number of threads")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    barcode_out = f"{args.output_prefix}.barcode_association.tsv"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)

    stats_out = f"{args.output_prefix}.barcode_association.stats.tsv"
    if os.path.exists(stats_out):
        os.remove(stats_out)

    print(f"Detecting variant and barcode associations, please wait...", flush=True)
    list_results = []
    for i, chunk_result in enumerate(process_long_reads_in_chunk(args.longread)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Processed chunk {i+1} with {args.chunk_size} read pairs", flush=True)
        if not chunk_result.is_empty():
            list_results.append(chunk_result)
    print(f"Finished processing all the reads.", flush=True)

    # -- clean and format the extracted barcodes from reads -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating barcode results, please wait ...", flush = True)
    list_results_filtered = [df for df in list_results if df.height > 0]
    if list_results_filtered:
        df_barcode = pl.concat(list_results_filtered, how = "vertical")
        df_barcode_counts = ( df_barcode.group_by(["barcode_seq", "vhh_seq", "target_seq"])
                                        .agg(pl.sum("count").alias("count")) )
    else:
        with open(barcode_out, "w") as f:
            f.write("no barcode found in the reads, please check your up/down-stream sequences!\n")
        exit(0)

    count_processed_reads = df_barcode_counts["count"].sum()
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Filtering barcode results, please wait ...", flush = True)

    # -- (1) barcode_seq == "upstream not found" -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode upstream match, please wait ...", flush = True)
    mask = df_barcode_counts["barcode_seq"] == "upstream not found"
    count_barup_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (2) barcode_seq == "downstream not found" -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode downstream match, please wait ...", flush = True)
    mask = df_barcode_counts["barcode_seq"] == "downstream not found"
    count_bardown_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (3) barcode pattern check -- #
    if args.barcode_temp is not None:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode template match, please wait ...", flush=True)
        df_barcode_counts = df_barcode_counts.with_columns(
            pl.col("barcode_seq").map_elements(
                lambda b: check_barcode(b, args.barcode_temp.upper(), args.barcode_mismatch),
                return_dtype=pl.Utf8
            ).alias("barcode_corrected")
        )
        count_barcode_template = int(df_barcode_counts.filter(pl.col("barcode_corrected").is_null())["count"].sum())
        df_barcode_counts = df_barcode_counts.filter(pl.col("barcode_corrected").is_not_null())
        df_barcode_counts = df_barcode_counts.drop("barcode_corrected")        

    # -- (4) vhh_seq == "upstream not found" -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against vhh upstream match, please wait ...", flush=True)
    mask = df_barcode_counts["vhh_seq"] == "upstream not found"
    count_vhhup_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (5) vhh_seq == "downstream not found" -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against vhh downstream match, please wait ...", flush=True)
    mask = df_barcode_counts["vhh_seq"] == "downstream not found"
    count_vhhdown_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (6) vhh length check -- #
    if args.vhh_len is not None:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against vhh length match, please wait ...", flush=True)
        vhh_lens = set(map(int, args.vhh_len.split(",")))
        mask = df_barcode_counts["vhh_seq"].str.len_chars().isin(vhh_lens)
        count_vhh_length = df_barcode_counts.filter(mask)["count"].sum()
        df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (7) target_seq == "upstream not found" -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against target upstream match, please wait ...", flush=True)
    mask = df_barcode_counts["target_seq"] == "upstream not found"
    count_tarup_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (8) target_seq == "downstream not found" -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against target downstream match, please wait ...", flush=True)
    mask = df_barcode_counts["target_seq"] == "downstream not found"
    count_tardown_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (9) target length check -- #
    if args.target_len is not None:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against target length match, please wait ...", flush=True)
        target_lens = set(map(int, args.target_len.split(",")))
        mask = df_barcode_counts["target_seq"].str.len_chars().isin(target_lens)
        count_target_length = df_barcode_counts.filter(mask)["count"].sum()
        df_barcode_counts = df_barcode_counts.filter(~mask)  
    
    # -- (10) count <= args.min_barcov -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode minimal coverage, please wait ...", flush=True)
    mask = df_barcode_counts["count"] < args.min_barcov
    count_low_barcov = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- (11) some barcodes match multiple targets, only keep the one with the highest count -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode multiple variant matches, please wait ...", flush=True)
    count_before_filter = df_barcode_counts["count"].sum()
    df_barcode_counts = ( df_barcode_counts.sort("count", descending = True)
                                           .group_by("barcode_seq")
                                           .agg([pl.first("vhh_seq").alias("vhh_seq"), 
                                                 pl.first("target_seq").alias("target_seq"), 
                                                 pl.first("count").alias("count")]) )

    # -- (12) remainning records -- #
    df_barcode_counts = df_barcode_counts.sort(["vhh_seq", "target_seq"])
    count_effective_reads = df_barcode_counts["count"].sum()

    # -- output results and stats -- #
    print(f"Creating output files, please wait...", flush=True)
    df_barcode_counts.write_csv(barcode_out, separator = "\t")

    with open(stats_out, "w") as f:
        f.write(f"Total reads processed: {count_processed_reads}\n")
        f.write(f"Total reads with barcode upstream not found: {count_barup_notfound}\n")
        f.write(f"Total reads with barcode downstream not found: {count_bardown_notfound}\n")
        if args.barcode_temp is not None:
            f.write(f"Total reads with barcode template mismatch: {count_barcode_template}\n")
        f.write(f"Total reads with vhh upstream not found: {count_vhhup_notfound}\n")
        f.write(f"Total reads with vhh downstream not found: {count_vhhdown_notfound}\n")
        if args.vhh_len is not None:
            f.write(f"Total reads with vhh length mismatch ({args.vhh_len}): {count_vhh_length}\n")
        f.write(f"Total reads with target upstream not found: {count_tarup_notfound}\n")
        f.write(f"Total reads with target downstream not found: {count_tardown_notfound}\n")
        if args.target_len is not None:
            f.write(f"Total reads with target length mismatch ({args.target_len}): {count_target_length}\n")
        f.write(f"Total reads with barcode coverage < {args.min_barcov}: {count_low_barcov}\n")
        f.write(f"Total reads with barcodes matching multiple variants: {count_before_filter - count_effective_reads}\n")
        f.write(f"Total effective reads: {count_effective_reads}\n")
