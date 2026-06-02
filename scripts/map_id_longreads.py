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
import edlib
import gzip
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice

#-------------------------------------------------------
# sequence processing functions
#-------------------------------------------------------
def detect_separator(filename):
    opener = gzip.open if filename.endswith(".gz") else open

    with opener(filename, "rt") as f:
        line = f.readline()

    tabs = line.count("\t")
    commas = line.count(",")

    if tabs > commas:
        return "\t"
    elif commas > tabs:
        return ","
    else:
        raise ValueError("Cannot determine delimiter")

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

#-------------------------------------------------------
# parallel processing functions
#-------------------------------------------------------
def process_record(record: tuple) -> list:
    """
    extract sequence by up/down markers in a record
    Parameters:
        -- record: tuple (barcode_seq, vhh_seq, target_seq) for record
    Returns:
        -- tuple: (barcode_seq, vhh_id, vhh_seq, target_id, target_seq)
    """
    barcode_seq = record[0]
    vhh_seq     = record[1]
    target_seq  = record[2]

    if vhh_seq == "" or target_seq == "":
        return (None, None, None, None, None)
    
    if len(vhh_seq) < vhh_min_len - 2 * args.max_mismatches or \
        len(vhh_seq) > vhh_max_len + 2 * args.max_mismatches or \
        len(target_seq) < target_min_len - 2 * args.max_mismatches or \
        len(target_seq) > target_max_len + 2 * args.max_mismatches:
        return (None, None, None, None, None)

    vhh_seq_rc    = reverse_complement(vhh_seq)
    target_seq_rc = reverse_complement(target_seq)

    vhh_id = None
    if vhh_seq in vhh_map:
        vhh_id = vhh_map[vhh_seq]
    elif vhh_seq_rc in vhh_map:
        vhh_id = vhh_map[vhh_seq_rc]
    else:
        for known_seq, known_id in vhh_map.items():
            found = match_hamming_numpy(vhh_seq, known_seq, args.max_mismatches)
            if found != -1:
                vhh_id = known_id
                break

            found = match_hamming_numpy(vhh_seq_rc, known_seq, args.max_mismatches)
            if found != -1:
                vhh_id = known_id
                break

            dist = edlib.align(vhh_seq, known_seq, mode = "NW", k = 2 * args.max_mismatches)["editDistance"]
            if dist != -1:
                vhh_id = known_id
                break

            dist = edlib.align(vhh_seq_rc, known_seq, mode = "NW", k = 2 * args.max_mismatches)["editDistance"]
            if dist != -1:
                vhh_id = known_id
                break
    
    target_id = None
    if target_seq in target_map:
        target_id = target_map[target_seq]
    elif target_seq_rc in target_map:
        target_id = target_map[target_seq_rc]
    else:
        for known_seq, known_id in target_map.items():
            found = match_hamming_numpy(target_seq, known_seq, args.max_mismatches)
            if found != -1:
                target_id = known_id
                break

            found = match_hamming_numpy(target_seq_rc, known_seq, args.max_mismatches)
            if found != -1:
                target_id = known_id
                break

            dist = edlib.align(target_seq, known_seq, mode = "NW", k = 2 * args.max_mismatches)["editDistance"]
            if dist != -1:
                target_id = known_id
                break

            dist = edlib.align(target_seq_rc, known_seq, mode = "NW", k = 2 * args.max_mismatches)["editDistance"]
            if dist != -1:
                target_id = known_id
                break

    return (barcode_seq, vhh_id, vhh_seq, target_id, target_seq)

def batch_process(batch: list) -> list:
    """
    Process a batch of long reads to extract sequences.
    Parameters:
        -- batch: list of records
    Returns:
        -- list of tuples
    """
    results = []
    for record in batch.iter_rows():
        result = process_record(record)
        results.append(result)
    return results

def function_processpool(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process(args)

def process_records_in_chunk(path_file: str):
    """
    Read the file in chunks and process records in parallel
    Parameters:
        -- path_file: path to the file
    Yields:
        -- DataFrame: records
    """
    reader = pl.read_csv_batched(
        path_file,
        separator = "\t",
        columns = ["barcode_seq", "vhh_seq", "target_seq"],
        has_header = True
    )
    
    with ProcessPoolExecutor(max_workers = args.threads) as executor:
        while True:
            record_chunk = reader.next_batches(args.chunk_size)
            if not record_chunk:
                break

            batch_size = min(args.chunk_size, 1000)
            record_batches = [
                pl.concat(record_chunk[i:i + batch_size])
                for i in range(0, len(record_chunk), batch_size)
            ]

            list_records = []
            futures = [ executor.submit(function_processpool, batch) for batch in record_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                list_records.append(pl.DataFrame(batch_result, schema = ["barcode_seq", "vhh_id", "vhh_seq", "target_id", "target_seq"], orient = "row"))

                # -- free memory -- #
                del batch_result
                gc.collect()

            df_yield = pl.concat(list_records, how = "vertical")

            # -- free memory -- #
            del record_chunk, record_batches, futures, list_records
            gc.collect()

            yield df_yield

#-------------------------------------------------------
# main function
#-------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Assign IDs to associated VNN and target sequences", allow_abbrev = False)
    parser.add_argument("--associate_file",   type = str, required = True,       help = "the path of long read FASTQ file")
    parser.add_argument("--vhh_id_file",      type = str, required = True,       help = "Upstream flank sequence of barcode in long read")
    parser.add_argument("--vhh_cols",         type = str, required = True,       help = "columns for vhh id and sequence in vhh_id_file, separated by comma, first is id")
    parser.add_argument("--target_id_file",   type = str, required = True,       help = "Downstream flank sequence of barcode in long read")
    parser.add_argument("--target_cols",      type = str, required = True,       help = "columns for target id and sequence in target_id_file, separated by comma, first is id")
    parser.add_argument("--max_mismatches",   type = int, default = 2,           help = "Max mismatches allowed in matches")
    parser.add_argument("--output_dir",       type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",    type = str, required = True,       help = "output prefix")
    parser.add_argument("--chunk_size",       type = int, default = 10000,       help = "Chunk size for processing reads")
    parser.add_argument("--threads",          type = int, default = 40,          help = "Number of threads")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    map_id_out = f"{args.output_prefix}.map_id.tsv"
    if os.path.exists(map_id_out):
        os.remove(map_id_out)

    # -- read inputs -- #
    print(f"Reading VHH file, please wait ...", flush = True)
    sep = detect_separator(args.vhh_id_file)
    cols = args.vhh_cols.split(",")
    vhh_ids = pl.read_csv(args.vhh_id_file, has_header = True, columns = cols, separator = sep)
    vhh_ids = vhh_ids.rename({ cols[0]: "vhh_id", cols[1]: "vhh_seq" })
    vhh_map = dict(zip(vhh_ids["vhh_seq"], vhh_ids["vhh_id"]))
    vhh_min_len = min(map(len, vhh_map))
    vhh_max_len = max(map(len, vhh_map))
    del vhh_ids
    gc.collect()

    print(f"Reading target file, please wait ...", flush = True)
    sep = detect_separator(args.target_id_file)
    cols = args.target_cols.split(",")
    target_ids = pl.read_csv(args.target_id_file, has_header = True, columns = cols, separator = sep)
    target_ids = target_ids.rename({ cols[0]: "target_id", cols[1]: "target_seq" })
    target_map = dict(zip(target_ids["target_seq"], target_ids["target_id"]))
    target_min_len = min(map(len, target_map))
    target_max_len = max(map(len, target_map))
    del target_ids
    gc.collect()

    print(f"Mapping IDs to records, please wait ...", flush = True)
    list_results = []
    for i, chunk_result in enumerate(process_records_in_chunk(args.associate_file)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Processed chunk {i+1} with {args.chunk_size} records", flush=True)
        if not chunk_result.is_empty():
            list_results.append(chunk_result)
    print(f"Finished processing all the records.", flush=True)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating barcode results, please wait ...", flush = True)
    list_results_filtered = [df for df in list_results if df.height > 0]
    if list_results_filtered:
        df_records = pl.concat(list_results_filtered, how = "vertical")
        df_records.write_csv(map_id_out, separator = "\t")
    else:
        with open(map_id_out, "w") as f:
            f.write("something wrong with input files, no records found!\n")
        exit(0)
