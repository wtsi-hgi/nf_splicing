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
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import islice

from sequence_utils import (
    pigz_open,
    base_cov_iter
)

#-- processing functions --#
def init_worker():
    """
    Initilize worker
    """
    global global_dict_exon_pos, global_dict_splicing_counts
    global_dict_exon_pos = dict_exon_pos
    global_dict_splicing_counts = dict_splicing_counts

def process_base_cov(base_cov: tuple) -> list:
    """
    Add canonical counts to base coverage
    Parameters:
        -- base_cov: tuple (var_id, base_pos, base_cov)
    Returns:
        -- tuple: tuple: (var_id, base_pos, base_cov) 
    """
    var_id, base_pos, base_cov = base_cov

    if var_id not in global_dict_exon_pos or \
       var_id not in global_dict_splicing_counts:
        return (None, None, None)

    e1_start, e1_end = global_dict_exon_pos[var_id][0]
    e2_start, e2_end = global_dict_exon_pos[var_id][1]
    e3_start, e3_end = global_dict_exon_pos[var_id][2]

    if base_pos < e1_start or base_pos > e3_end:
        return (None, None, None)

    inclusion_count, skipping_count = global_dict_splicing_counts[var_id]

    if e1_start <= base_pos <= e1_end:
        base_cov += inclusion_count + skipping_count
    elif e2_start <= base_pos <= e2_end:
        base_cov += inclusion_count
    elif e3_start <= base_pos <= e3_end:
        base_cov += inclusion_count + skipping_count

    return (var_id, base_pos, base_cov)

def batch_process(batch_base_covs: list) -> list:
    """
    Process a batch of base coverages to extract variants and barcodes.
    Parameters:
        -- batch_base_covs
    Returns:
        -- list of tuples
    """
    results = []
    for base_cov in batch_base_covs:
        result = process_base_cov(base_cov)
        results.append(result)
    return results

def function_processpool(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process(args)

def update_base_cov_in_chunks(base_cov_file: str, chunk_size: int, threads: int):
    """
    Read base coverage in chunks and adding counts
    Parameters:
        -- base_cov_file: base coverage file (for novel splicing events)
        -- chunk_size: Number of base_covs to process in each chunk
        -- threads: Number of threads to use for processing
    Returns:
        Generator yielding updated base coverage and splicing site usage (SSU)
    """
    fh_base_cov = io.TextIOWrapper(pigz_open(base_cov_file).stdout) if base_cov_file.endswith(".gz") else open(base_cov_file)
    base_cov_tuple = base_cov_iter(fh_base_cov)

    with ProcessPoolExecutor(max_workers = args.threads, initializer = init_worker) as executor:
        while True:
            base_cov_chunk = list(islice(base_cov_tuple, args.chunk_size))
            if not base_cov_chunk:
                break

            batch_size = max(2000, args.chunk_size // (args.threads * 4))
            base_cov_batches = [
                base_cov_chunk[i:i+batch_size]
                for i in range(0, len(base_cov_chunk), batch_size)
            ]

            list_base_covs = []
            futures = [ executor.submit(function_processpool, batch) for batch in base_cov_batches ]
            for future in as_completed(futures):
                batch_result = future.result()
                list_base_covs.append(pl.DataFrame(batch_result, schema = ["var_id", "base_pos", "base_cov"], orient = "row"))

            df_yields = pl.concat(list_base_covs, how = "vertical").filter(pl.col("var_id").is_not_null())
            yield df_yields
    fh_base_cov.close()            

#-- duckdb functions --#
def duckdb_merge(chunk_files: list, tmp_dir: str) -> pl.DataFrame:
    """
    Merge all chunk parquet files using DuckDB with automatic spill-to-disk.
    Parameters:
        -- chunk_files: list of parquet file paths to merge
        -- tmp_dir:     directory to write the final merged parquet and spill files
    Returns:
        -- pl.DataFrame with columns ["var_id", "base_pos", "base_cov", "max_cov", "base_ssu"]
    """
    con = duckdb.connect()
    con.execute(f"SET temp_directory='{tmp_dir}'")
    con.execute(f"SET memory_limit='{args.db_mem_limit}'")
    con.execute(f"SET threads={args.threads}")

    file_list = ", ".join(f"'{f}'" for f in chunk_files)
    final_path = os.path.join(tmp_dir, "final.parquet")

    con.execute(f"""
        COPY (
            WITH coverage AS (
                SELECT 
                    var_id, 
                    base_pos, 
                    base_cov, 
                    MAX(base_cov) OVER(
                        PARTITION BY var_id
                    ) AS max_cov
                FROM read_parquet([{file_list}]) 
            )
            SELECT 
                var_id, 
                base_pos, 
                base_cov, 
                max_cov, 
                CASE
                    WHEN max_cov = 0 THEN 0
                    ELSE ROUND(base_cov / max_cov, 2) 
                END AS base_ssu
            FROM coverage
            ORDER BY var_id, base_pos
        )
        TO '{final_path}' (FORMAT PARQUET)
    """)
    con.close()

    return pl.read_parquet(final_path)

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Extract variant and barcode from paired-end FASTQ files.", allow_abbrev = False)
    parser.add_argument("--exon_pos",         type = str,       required = True,       help = "Exon position file")
    parser.add_argument("--splicing_counts",  type = str,       required = True,       help = "Splicing counts file (for canonical splicing events)")
    parser.add_argument("--base_cov",         type = str,       required = True,       help = "Base coverage file (for novel splicing events)")
    parser.add_argument("--output_dir",       type = str,       default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",    type = str,       required = True,       help = "output prefix")
    parser.add_argument("--resume_tmp",       action = "store_true",                   help = "Whether to resume the process and keep temporary files")
    parser.add_argument("--chunk_size",       type = int,       default = 100000,      help = "Chunk size for processing reads")
    parser.add_argument("--threads",          type = int,       default = 20,          help = "Number of threads")
    parser.add_argument("--db_mem_limit",     type = str,       default = "10GB",      help = "Memory limit for DuckDB during merging")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    #-- input file --#
    print(f"Reading input files, please wait...", flush=True)
    exon_pos = pl.read_csv(args.exon_pos, separator = "\t", has_header = False, new_columns = ["var_id", "exon_id", "exon_start", "exon_end"])
    splicing_counts = pl.read_csv(args.splicing_counts, separator = "\t", has_header = True, columns = ["var_id", "canonical_inclusion", "canonical_skipping"])

    dict_exon_pos = {}
    for row in exon_pos.iter_rows(named = True):
        var_id = row["var_id"]

        # inintialise
        if var_id not in dict_exon_pos:
            dict_exon_pos[var_id] = [None, None, None]

        match row["exon_id"]:
            case "E1":
                dict_exon_pos[var_id][0] = (row["exon_start"], row["exon_end"])
            case "E2":
                dict_exon_pos[var_id][1] = (row["exon_start"], row["exon_end"])
            case "E3":
                dict_exon_pos[var_id][2] = (row["exon_start"], row["exon_end"])
            case _:
                pass

    dict_splicing_counts = { row["var_id"]: (row["canonical_inclusion"], row["canonical_skipping"]) for row in splicing_counts.iter_rows(named = True) }

    del exon_pos
    del splicing_counts
    gc.collect()

    #-- output file --#
    ssu_out = f"{args.output_prefix}.splicing_ssu.tsv"
    if os.path.exists(ssu_out):
        os.remove(ssu_out)

    #-- processing --#
    print(f"Processing, please wait...", flush=True)

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
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} No existing chunk files found in {tmp_dir}, starting from scratch.", flush=True)
        else:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Found {len(chunk_files)} existing chunk files in {tmp_dir}, resuming from these files.", flush=True)

    if not chunk_files:
        for i, chunk_result in enumerate(update_base_cov_in_chunks(args.base_cov, args.chunk_size, args.threads)):
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> Processed chunk {i+1} with {args.chunk_size} lines", flush=True)
            if not chunk_result.is_empty():
                tmp_path = os.path.join(tmp_dir, f"tmp_chunk_{i}.parquet")
                chunk_result.write_parquet(tmp_path)
                chunk_files.append(tmp_path)    

    if not chunk_files:
        with open(ssu_out, "w") as f:
            f.write("no data found in the base coverage file, please check your inputs!\n")
        exit(0)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating SSU results, please wait...", flush=True)
    df_base_cov = ( 
        duckdb_merge(chunk_files, tmp_dir)
        .with_columns(
            pl.col(["base_pos", "base_cov", "max_cov"]).cast(pl.Int64),
            pl.col("base_ssu").cast(pl.Float64)
        )
    )

    if not args.resume_tmp:
        for f in chunk_files:
            os.remove(f)
        os.remove(os.path.join(tmp_dir, "final.parquet"))
        os.rmdir(tmp_dir)

    df_base_cov.write_csv(ssu_out, separator = "\t")    
