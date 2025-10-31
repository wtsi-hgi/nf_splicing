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
    fastq_iter_pe,
    extract_sequence,
    reverse_complement,
    check_barcode
)

#-- functions --#
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

    variant_seq = extract_sequence(read1_seq, args.variant_up, args.variant_down, args.max_mismatches)
    barcode_seq = extract_sequence(read2_seq, args.barcode_up, args.barcode_down, args.max_mismatches)
    
    return (variant_seq, barcode_seq)

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

def function_processpool_pe(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_pe_pairs(args)
    
def process_pe_pairs_in_chunk(path_read1, path_read2):
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

    with ProcessPoolExecutor(max_workers = args.threads) as executor:
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

            list_barcodes = []
            futures = [ executor.submit(function_processpool_pe, batch) for batch in read_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                list_barcodes.append(pl.DataFrame(batch_result, schema=["variant_seq", "barcode_seq"]))

                # -- free memory -- #
                del batch_result
                gc.collect()

            df_yield = ( pl.concat(list_barcodes, how="vertical")
                           .group_by(["variant_seq", "barcode_seq"])
                           .agg(pl.count()) )

            # -- free memory -- #
            del read_chunk, read_batches, futures, list_barcodes
            gc.collect()

            yield df_yield
    fh_read1.close()
    fh_read2.close()

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Extract variant and barcode from paired-end FASTQ files.", allow_abbrev = False)
    parser.add_argument("--read1",            type = str, required = True,       help = "Read 1 FASTQ file")
    parser.add_argument("--read2",            type = str, required = True,       help = "Read 2 FASTQ file")
    parser.add_argument("--variant_up",       type = str, required = True,       help = "Upstream flank sequence of variant in read1")
    parser.add_argument("--variant_down",     type = str, required = True,       help = "Downstream flank sequence of variant in read1")
    parser.add_argument("--variant_check",    action="store_true",               help = "Enable variant checking against length")
    parser.add_argument("--variant_len",      type = int,                        help = "Length of the variant sequence")
    parser.add_argument("--barcode_up",       type = str, required = True,       help = "Upstream flank sequence of barcode in read2")
    parser.add_argument("--barcode_down",     type = str, required = True,       help = "Downstream flank sequence of barcode in read2")
    parser.add_argument("--barcode_check",    action="store_true",               help = "Enable barcode checking against template")
    parser.add_argument("--barcode_temp",     type = str,                        help = "Template for barcode sequence (e.g., 'NNATNNNNATNNNNATNNNN')")
    parser.add_argument("--barcode_mismatch", type = int, default = 1,           help = "Number of mismatches allowed in barcode checking")
    parser.add_argument("--max_mismatches",   type = int, default = 2,           help = "Max mismatches allowed in up/down matches")
    parser.add_argument("--min_barcov",       type = int, default = 2,           help = "Minimum coverage for barcode-variant association")
    parser.add_argument("--output_dir",       type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",    type = str, required = True,       help = "output prefix")
    parser.add_argument("--chunk_size",       type = int, default = 100000,      help = "Chunk size for processing reads")
    parser.add_argument("--threads",          type = int, default = 40,          help = "Number of threads")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.variant_check and args.variant_len is None:
        parser.error("--variant_len is required when --variant_check is set")

    if args.barcode_check and args.barcode_temp is None:
        parser.error("--barcode_temp is required when --barcode_check is set")

    barcode_out = f"{args.output_prefix}.barcode_association.tsv"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)

    stats_out = f"{args.output_prefix}.barcode_association.stats.tsv"
    if os.path.exists(stats_out):
        os.remove(stats_out)

    print(f"Detecting variant and barcode associations, please wait...", flush=True)
    list_results = []
    for i, chunk_result in enumerate(process_pe_pairs_in_chunk(args.read1, args.read2)):
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Processed chunk {i+1} with {args.chunk_size} read pairs", flush=True)
        if not chunk_result.is_empty():
            list_results.append(chunk_result)
    print(f"Finished processing all the reads.", flush=True)

    # -- clean and format the extracted barcodes from reads -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating barcode results, please wait ...", flush = True)
    list_results_filtered = [df for df in list_results if df.height > 0]
    if list_results_filtered:
        df_barcode = pl.concat(list_results_filtered, how = "vertical")
        df_barcode_counts = ( df_barcode.group_by(["variant_seq", "barcode_seq"])
                                        .agg(pl.sum("count").alias("count")) )
    else:
        with open(barcode_out, "w") as f:
            f.write("no barcode found in the reads, please check your barcode marker or template!\n")
        exit(0)

    count_processed_reads = df_barcode_counts["count"].sum()
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Filtering barcode results, please wait ...", flush = True)

    # --- (1) variant_seq == "upstream not found"
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against variant upstream match, please wait ...", flush=True)
    mask = df_barcode_counts["variant_seq"] == "upstream not found"
    count_varup_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # --- (2) variant_seq == "downstream not found"
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against variant downstream match, please wait ...", flush=True)
    mask = df_barcode_counts["variant_seq"] == "downstream not found"
    count_vardown_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # --- (3) barcode_seq == "upstream not found"
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode upstream match, please wait ...", flush=True)
    mask = df_barcode_counts["barcode_seq"] == "upstream not found"
    count_barup_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # --- (4) barcode_seq == "downstream not found"
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode downstream match, please wait ...", flush=True)
    mask = df_barcode_counts["barcode_seq"] == "downstream not found"
    count_bardown_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # --- (5) len(variant_seq) < args.variant_len
    if args.variant_check:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against variant length match, please wait ...", flush=True)
        mask = df_barcode_counts["variant_seq"].str.len_chars() < args.variant_len
        count_variant_length = df_barcode_counts.filter(mask)["count"].sum()
        df_barcode_counts = df_barcode_counts.filter(~mask)

    # --- (6) barcode pattern check
    if args.barcode_check:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode template match, please wait ...", flush=True)
        rows = []
        for variant_seq, barcode_seq, count in df_barcode_counts.iter_rows():
            barcode_corrected = check_barcode(barcode_seq, args.barcode_temp.upper(), args.barcode_mismatch)
            if barcode_corrected is None:
                count_barcode_template += count
            else:
                rows.append((variant_seq, barcode_seq, count))
        df_barcode_counts = pl.DataFrame(rows, schema=["variant_seq", "barcode_seq", "count"])

    # --- (7) count <= args.min_barcov
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode minimal coverage, please wait ...", flush=True)
    mask = df_barcode_counts["count"] <= args.min_barcov
    count_low_barcov = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # --- (8) remaining records
    count_effective_reads = df_barcode_counts["count"].sum()

    print(f"Creating output files, please wait...", flush=True)
    df_barcode_counts.write_csv(barcode_out, separator = "\t")

    with open(stats_out, "w") as f:
        f.write(f"Total reads processed: {len(count_processed_reads)}\n")
        f.write(f"Total reads with variant upstream not found: {count_varup_notfound}\n")
        f.write(f"Total reads with variant downstream not found: {count_vardown_notfound}\n")
        f.write(f"Total reads with barcode upstream not found: {count_barup_notfound}\n")
        f.write(f"Total reads with barcode downstream not found: {count_bardown_notfound}\n")
        if args.variant_check:
            f.write(f"Total reads with variant length inconsistent: {count_variant_length}\n")
        if args.barcode_check:
            f.write(f"Total reads with barcode pattern mismatch: {count_barcode_template}\n")
        f.write(f"Total effective reads: {count_effective_reads}\n")
    

