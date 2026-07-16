#-- import modules --#
import io
import os
import sys
import argparse
import re
import gc
import subprocess
import duckdb
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice

from sequence_utils import (
    pigz_open,
    fastq_iter_pe,
    extract_sequence,
    reverse_complement,
    check_barcode
)

#-------------------------------------------------------
# parallel processing functions for paired-end reads
#-------------------------------------------------------
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

    target_seq = read1_seq
    if args.barcode_up is None and args.barcode_down is None:
        barcode_seq = read2_seq
    else:
        barcode_seq = extract_sequence(read2_seq, args.barcode_up, args.barcode_down, args.max_mismatches)
    
    return (target_seq, barcode_seq)

def batch_process_pe_pairs(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract targets and barcodes.
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
            batch_size = min(args.chunk_size, 40000)
            read_batches = [
                read_chunk[i:i+batch_size]
                for i in range(0, len(read_chunk), batch_size)
            ]

            list_barcodes = []
            futures = [ executor.submit(function_processpool_pe, batch) for batch in read_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                list_barcodes.append(pl.DataFrame(batch_result, schema = ["target_seq", "barcode_seq"], orient = "row"))

                # -- free memory -- #
                del batch_result
                gc.collect()

            df_yield = ( pl.concat(list_barcodes, how = "vertical")
                           .group_by(["target_seq", "barcode_seq"])
                           .agg(pl.len().alias("count")) )

            # -- free memory -- #
            del read_chunk, read_batches, futures, list_barcodes
            gc.collect()

            yield df_yield
    fh_read1.close()
    fh_read2.close()

#-------------------------------------------------------
# memory-efficient data merge
#-------------------------------------------------------
def duckdb_merge(chunk_files: list, tmp_dir: str) -> pl.DataFrame:
    """
    Merge all chunk parquet files using DuckDB with automatic spill-to-disk.
    Parameters:
        -- chunk_files: list of parquet file paths to merge
        -- tmp_dir:     directory to write the final merged parquet and spill files
    Returns:
        -- pl.DataFrame with columns ["target_seq", "barcode_seq", "count"]
    """
    con = duckdb.connect()
    con.execute(f"SET temp_directory='{tmp_dir}'")
    con.execute(f"SET memory_limit='{args.db_mem_limit}'")
    con.execute(f"SET threads={args.threads}") 

    file_list = ", ".join(f"'{f}'" for f in chunk_files)
    final_path = os.path.join(tmp_dir, "final.parquet")

    con.execute(f"""
        COPY (
            SELECT target_seq, barcode_seq, SUM(count) AS count
            FROM read_parquet([{file_list}])
            GROUP BY target_seq, barcode_seq
            ORDER BY target_seq
        )
        TO '{final_path}' (FORMAT PARQUET)
    """)
    con.close()

    return pl.read_parquet(final_path)

#-------------------------------------------------------
# create histogram of barcode counts
#-------------------------------------------------------
def create_barcode_count_histogram(df_barcode_counts: pl.DataFrame, output_path: str):
    """
    Create a histogram of barcode counts and save to file.
    Parameters:
        -- df_barcode_counts: DataFrame with columns ["target_seq", "barcode_seq", "count"]
        -- output_path: path to save the histogram plot
    """
    counts = df_barcode_counts["count"].to_numpy()

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(counts, bins = 60, color = "royalblue", edgecolor = "white", linewidth = 0.4)
    ax.set_title("Histogram of Barcode Counts", fontsize = 14)
    ax.set_xlabel("Count", fontsize = 12)
    ax.set_ylabel("Number of barcodes", fontsize = 12)

    median_val = np.median(counts)
    mean_val = np.mean(counts)
    n_barcodes = len(counts)
    ax.axvline(median_val, color = "red", linestyle = "--", linewidth = 1.2, label = f"Median = {median_val:.1f}")
    ax.axvline(mean_val, color = "orange", linestyle = "--", linewidth = 1.2, label = f"Mean = {mean_val:.1f}")
    ax.legend(fontsize = 10)

    ax.text(0.98, 0.97, f"n barcodes = {n_barcodes:,}", transform = ax.transAxes, ha = "right", va = "top", fontsize = 8, color = "black")
 
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    fig.savefig(output_path, dpi = 150)
    plt.close(fig)

#-------------------------------------------------------
# main execution
#-------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Extract target and barcode from paired-end FASTQ files.", allow_abbrev = False)
    parser.add_argument("--read1",             type = str, required = True,       help = "Read 1 FASTQ file")
    parser.add_argument("--read2",             type = str, required = True,       help = "Read 2 FASTQ file")
    parser.add_argument("--barcode_up",        type = str, default = None,        help = "Upstream flank sequence of barcode in read2")
    parser.add_argument("--barcode_down",      type = str, default = None,        help = "Downstream flank sequence of barcode in read2")
    parser.add_argument("--max_mismatches",    type = int, default = 2,           help = "Max mismatches allowed in up/down matches")
    parser.add_argument("--barcode_len",       type = str, default = None,        help = "Length of the barcode sequence")
    parser.add_argument("--restrict_site",     type = str, default = None,        help = "Sequence of restriction enzyme cutting site in the vector")
    parser.add_argument("--restrict_mismatch", type = int, default = 1,           help = "Number of mismatches allowed in restriction site checking")
    parser.add_argument("--min_barcov",        type = int, default = 2,           help = "Minimum coverage for barcode-target association")
    parser.add_argument("--resume_tmp",        action = "store_true",             help = "Whether to resume the process and keep temporary files")
    parser.add_argument("--output_dir",        type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",     type = str, required = True,       help = "output prefix")
    parser.add_argument("--chunk_size",        type = int, default = 100000,      help = "Chunk size for processing reads")
    parser.add_argument("--threads",           type = int, default = 40,          help = "Number of threads")
    parser.add_argument("--db_mem_limit",      type = str, default = "60GB",      help = "Memory limit for DuckDB during merging")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    # -- creating outputs -- #
    barcode_out = f"{args.output_prefix}.barcode_association.tsv"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)

    stats_out = f"{args.output_prefix}.barcode_association.stats.tsv"
    if os.path.exists(stats_out):
        os.remove(stats_out)

    #-- processing --#
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Detecting target and barcode associations, please wait...", flush=True)

    tmp_dir = os.path.join(args.output_dir, args.output_prefix + "_tmp")
    if os.path.exists(tmp_dir):
        if not args.resume_tmp:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Warning: Temporary directory {tmp_dir} already exists, it will be removed and recreated.", flush=True)
            for f in os.listdir(tmp_dir):
                os.remove(os.path.join(tmp_dir, f))
        else:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Resuming from existing temporary directory {tmp_dir}.", flush=True)
    else:
        os.makedirs(tmp_dir, exist_ok = True)
    
    chunk_files = []
    if args.resume_tmp:
        existing = {f for f in os.listdir(tmp_dir) if f.startswith("tmp_chunk_") and f.endswith(".parquet")}
        chunk_files = [os.path.join(tmp_dir, f) for f in sorted(existing)]
        if not chunk_files:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} No existing chunk files found in {tmp_dir}, starting fresh.", flush=True)
        else:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Found {len(chunk_files)} existing chunk files in {tmp_dir}, resuming from these files.", flush=True)
    
    if not chunk_files:
        for i, chunk_result in enumerate(process_pe_pairs_in_chunk(args.read1, args.read2)):
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Processed chunk {i+1} with {args.chunk_size} read pairs", flush=True)
            if not chunk_result.is_empty():
                tmp_path = os.path.join(tmp_dir, f"tmp_chunk_{i}.parquet")
                chunk_result.write_parquet(tmp_path)
                chunk_files.append(tmp_path)
            del chunk_result
            gc.collect()

    if not chunk_files:
        with open(barcode_out, "w") as f:
            f.write("no barcode found in the reads, please check your barcode marker or template!\n")
        exit(0)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating barcode results, please wait...", flush=True)
    df_barcode_counts = duckdb_merge(chunk_files, tmp_dir).with_columns(pl.col("count").cast(pl.Int64))

    if not args.resume_tmp:
        for f in chunk_files:
            os.remove(f)
        os.rmdir(tmp_dir)

    count_processed_reads = df_barcode_counts["count"].sum()

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Filtering barcode results, please wait ...", flush = True)

    # -- barcode upstream not found -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode upstream match, please wait ...", flush=True)
    mask = df_barcode_counts["barcode_seq"] == "upstream not found"
    count_barup_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- barcode downstream not found -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode downstream match, please wait ...", flush=True)
    mask = df_barcode_counts["barcode_seq"] == "downstream not found"
    count_bardown_notfound = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- barcode length check -- #
    if args.barcode_len is not None:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode length, please wait ...", flush=True)
        valid_lengths = set(int(x.strip()) for x in str(args.barcode_len).split(","))
        mask = ~df_barcode_counts["barcode_seq"].str.len_chars().is_in(valid_lengths)
        count_barcode_length = df_barcode_counts.filter(mask)["count"].sum()
        df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- restriction site check -- #
    if args.restrict_site is not None:
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against restriction site match, please wait ...", flush=True)
        mask = df_barcode_counts["barcode_seq"].map_elements(
            lambda b: check_barcode(b, args.restrict_site.upper(), args.restrict_mismatch),
            return_dtype=pl.Boolean
        )
        count_barcode_restrict = df_barcode_counts.filter(~mask)["count"].sum()
        df_barcode_counts = df_barcode_counts.filter(mask)
    
    # -- barcode coverage check -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode minimal coverage, please wait ...", flush=True)
    mask = df_barcode_counts["count"] < args.min_barcov
    count_low_barcov = df_barcode_counts.filter(mask)["count"].sum()
    df_barcode_counts = df_barcode_counts.filter(~mask)

    # -- some barcodes match multiple targets, only keep the one with the highest count -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Filtering against barcode multiple target matches, please wait ...", flush=True)
    count_before_filter = df_barcode_counts["count"].sum()
    df_barcode_counts = ( df_barcode_counts.sort("count", descending = True)
                                           .group_by("barcode_seq")
                                           .agg([pl.first("target_seq").alias("target_seq"), pl.first("count").alias("count")]) )

    # -- remaining records -- #
    df_barcode_counts = df_barcode_counts.sort("target_seq")
    count_effective_reads = df_barcode_counts["count"].sum()

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Creating output files, please wait...", flush=True)
    df_barcode_counts.write_csv(barcode_out, separator = "\t")

    # -- statistics -- #
    n_effective_targets = df_barcode_counts["target_seq"].n_unique()
    n_avg_barcode_per_target = ( 
        df_barcode_counts.group_by("target_seq")
                         .agg(pl.count("barcode_seq").alias("barcode_count"))
                         .select(pl.col("barcode_count").mean())
                         .item() 
    )

    # -- write statistics -- #
    with open(stats_out, "w") as f:
        f.write(f"Total reads processed: {count_processed_reads}\n")
        f.write(f"Total reads with barcode upstream not found: {count_barup_notfound}\n")
        f.write(f"Total reads with barcode downstream not found: {count_bardown_notfound}\n")
        if args.barcode_len is not None:
            f.write(f"Total reads with barcode length inconsistent: {count_barcode_length}\n")
        if args.restrict_site is not None:
            f.write(f"Total reads with restriction site mismatch: {count_barcode_restrict}\n")
        f.write(f"Total reads with barcode coverage < {args.min_barcov}: {count_low_barcov}\n")
        f.write(f"Total reads with barcodes matching multiple targets: {count_before_filter - count_effective_reads}\n")
        f.write(f"Total effective reads: {count_effective_reads}\n")
        f.write(f"Total effective targets: {n_effective_targets}\n")
        f.write(f"Average number of barcodes per target: {n_avg_barcode_per_target:.2f}\n")

    # -- create histogram of barcode counts -- #
    hist_out = f"{args.output_prefix}.barcode_count_histogram.png"
    create_barcode_count_histogram(df_barcode_counts, hist_out)
