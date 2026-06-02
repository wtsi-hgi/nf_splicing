#-- import modules --#
import io
import os
import sys
import argparse
import gc
import gzip
import numpy as np
import polars as pl
import edlib
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

#-------------------------------------------------------
# sequence processing functions
#-------------------------------------------------------
def detect_separator(filename):
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(filename, "rt") as f:
        line = f.readline()
    tabs   = line.count("\t")
    commas = line.count(",")
    if tabs > commas:
        return "\t"
    elif commas > tabs:
        return ","
    else:
        raise ValueError("Cannot determine delimiter")

def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

# ---------------------------------------------------------------------------
# Vectorised matcher: compare ONE query against ALL known seqs at once
# ---------------------------------------------------------------------------
def best_hamming_match_vectorised(query: str,
                                   known_arrays: list,   # list of np.uint8 arrays, all same length k
                                   known_ids: list,
                                   max_mismatches: int):
    """
    Compare query against every known sequence in one numpy call.
    Returns the id of the first match within max_mismatches, or "".
    All sequences in known_arrays must be exactly len(query) long.
    """
    k = len(query)
    if not known_arrays:
        return ""
    q_arr = np.frombuffer(query.encode("ascii"), dtype=np.uint8)
    # Stack all known seqs → shape (N, k)
    mat = np.stack(known_arrays)                 # (N, k)
    mismatches = np.sum(mat != q_arr, axis=1)    # (N,)
    hits = np.where(mismatches <= max_mismatches)[0]
    if hits.size > 0:
        return known_ids[hits[0]]
    return ""

def best_edlib_match(query: str, known_seqs: list, known_ids: list, max_ed: int):
    """
    Try edlib (edit distance) against all known seqs, return first hit id or "".
    Only called when hamming already failed.
    """
    for seq, sid in zip(known_seqs, known_ids):
        dist = edlib.align(query, seq, mode="NW", k=max_ed)["editDistance"]
        if dist != -1:
            return sid
    return ""

# ---------------------------------------------------------------------------
# Worker initialiser – runs ONCE per worker process (no repeated pickling)
# ---------------------------------------------------------------------------
_vhh_map    = None
_target_map = None
_vhh_arrays = None   # list of np.uint8 per known vhh sequence
_vhh_ids_list   = None
_vhh_seqs_list  = None
_target_arrays  = None
_target_ids_list  = None
_target_seqs_list = None
_vhh_min_len    = None
_vhh_max_len    = None
_target_min_len = None
_target_max_len = None
_max_mismatches = None

def _worker_init(vhh_map, target_map, max_mismatches):
    """Called once when each worker process starts."""
    global _vhh_map, _target_map
    global _vhh_arrays, _vhh_ids_list, _vhh_seqs_list
    global _target_arrays, _target_ids_list, _target_seqs_list
    global _vhh_min_len, _vhh_max_len, _target_min_len, _target_max_len
    global _max_mismatches

    _vhh_map    = vhh_map
    _target_map = target_map
    _max_mismatches = max_mismatches

    # Pre-build numpy arrays for hamming vectorisation
    _vhh_seqs_list = list(vhh_map.keys())
    _vhh_ids_list  = list(vhh_map.values())
    _vhh_arrays    = [np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in _vhh_seqs_list]

    _target_seqs_list = list(target_map.keys())
    _target_ids_list  = list(target_map.values())
    _target_arrays    = [np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in _target_seqs_list]

    _vhh_min_len    = min(len(s) for s in _vhh_seqs_list)
    _vhh_max_len    = max(len(s) for s in _vhh_seqs_list)
    _target_min_len = min(len(s) for s in _target_seqs_list)
    _target_max_len = max(len(s) for s in _target_seqs_list)

# ---------------------------------------------------------------------------
# Per-record matching logic
# ---------------------------------------------------------------------------
def _resolve_id(query_seq: str,
                exact_map: dict,
                known_arrays: list,
                known_ids_list: list,
                known_seqs_list: list,
                max_mismatches: int) -> str:
    """
    Try exact → RC exact → vectorised hamming → edlib for both orientations.
    Returns matched id or "".
    """
    # 1. exact forward
    if query_seq in exact_map:
        return exact_map[query_seq]

    # 2. exact reverse complement
    rc = reverse_complement(query_seq)
    if rc in exact_map:
        return exact_map[rc]

    # 3. Hamming only makes sense when lengths match – filter by length
    q_len = len(query_seq)
    fwd_arrays = [a for a, s in zip(known_arrays, known_seqs_list) if len(s) == q_len]
    fwd_ids    = [i for i, s in zip(known_ids_list, known_seqs_list) if len(s) == q_len]

    rc_arrays  = fwd_arrays   # same length filter applies to rc (same len)
    rc_ids     = fwd_ids

    # vectorised hamming – forward
    hit = best_hamming_match_vectorised(query_seq, fwd_arrays, fwd_ids, max_mismatches)
    if hit:
        return hit

    # vectorised hamming – RC
    hit = best_hamming_match_vectorised(rc, rc_arrays, rc_ids, max_mismatches)
    if hit:
        return hit

    # 4. edlib (only if hamming found nothing)
    hit = best_edlib_match(query_seq, known_seqs_list, known_ids_list, 2 * max_mismatches)
    if hit:
        return hit

    hit = best_edlib_match(rc, known_seqs_list, known_ids_list, 2 * max_mismatches)
    return hit  # "" if nothing found

def process_record(record: tuple):
    barcode_seq, vhh_seq, target_seq = record

    if not vhh_seq or not target_seq:
        return None

    # Length pre-filter (fast reject before any alignment)
    mm = _max_mismatches
    if (len(vhh_seq) < _vhh_min_len - mm or
            len(vhh_seq) > _vhh_max_len + mm or
            len(target_seq) < _target_min_len - mm or
            len(target_seq) > _target_max_len + mm):
        return None

    vhh_id = _resolve_id(
        vhh_seq, _vhh_map, _vhh_arrays, _vhh_ids_list, _vhh_seqs_list, mm
    )
    target_id = _resolve_id(
        target_seq, _target_map, _target_arrays, _target_ids_list, _target_seqs_list, mm
    )

    return (barcode_seq, vhh_id, vhh_seq, target_id, target_seq)

# ---------------------------------------------------------------------------
# Batch processor (runs inside each worker)
# ---------------------------------------------------------------------------
def batch_process(rows: list) -> list:
    """rows is a plain Python list of tuples – no polars overhead in workers."""
    results = []
    for row in rows:
        result = process_record(row)
        if result is not None:
            results.append(result)
    return results

# ---------------------------------------------------------------------------
# Main reader / dispatcher
# ---------------------------------------------------------------------------
def process_records_in_chunk(path_file: str, args, vhh_map: dict, target_map: dict):
    """
    Read the TSV in streaming batches and dispatch to worker pool.
    Workers are initialised ONCE via _worker_init so maps are never re-pickled.
    """
    SCHEMA = ["barcode_seq", "vhh_id", "vhh_seq", "target_id", "target_seq"]

    reader = pl.read_csv_batched(
        path_file,
        separator="\t",
        columns=["barcode_seq", "vhh_seq", "target_seq"],
        has_header=True,
    )

    # Each worker gets maps once at startup
    initargs = (vhh_map, target_map, args.max_mismatches)

    with ProcessPoolExecutor(
        max_workers=args.threads,
        initializer=_worker_init,
        initargs=initargs,
    ) as executor:
        while True:
            # Pull a larger number of polars batches at once to keep workers fed
            record_chunk = reader.next_batches(args.threads * 2)
            if not record_chunk:
                break

            # Flatten into one DataFrame then split into per-worker row lists
            combined = pl.concat(record_chunk, how="vertical")
            rows_all = combined.rows()   # list of plain Python tuples – cheap to pickle

            # Split into equal sub-lists, one per worker thread
            chunk_size = max(1, len(rows_all) // args.threads)
            batches = [
                rows_all[i : i + chunk_size]
                for i in range(0, len(rows_all), chunk_size)
            ]

            futures = [executor.submit(batch_process, b) for b in batches]
            list_records = []
            for future in as_completed(futures):
                batch_result = future.result()
                if batch_result:
                    list_records.append(
                        pl.DataFrame(batch_result, schema=SCHEMA, orient="row")
                    )
                del batch_result

            if list_records:
                df_yield = pl.concat(list_records, how="vertical")
                df_yield = df_yield.filter(pl.any_horizontal(pl.all().is_not_null()))
            else:
                df_yield = pl.DataFrame([], schema=SCHEMA)

            del record_chunk, combined, rows_all, batches, futures, list_records
            gc.collect()

            yield df_yield

#-------------------------------------------------------
# main
#-------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assign IDs to associated VHH and target sequences",
        allow_abbrev=False,
    )
    parser.add_argument("--associate_file",  type=str, required=True)
    parser.add_argument("--vhh_id_file",     type=str, required=True)
    parser.add_argument("--vhh_cols",        type=str, required=True,
                        help="comma-separated: id_col,seq_col")
    parser.add_argument("--target_id_file",  type=str, required=True)
    parser.add_argument("--target_cols",     type=str, required=True,
                        help="comma-separated: id_col,seq_col")
    parser.add_argument("--max_mismatches",  type=int, default=2)
    parser.add_argument("--output_dir",      type=str, default=os.getcwd())
    parser.add_argument("--output_prefix",   type=str, required=True)
    parser.add_argument("--threads",         type=int, default=40)

    args, unknown = parser.parse_known_args()
    if unknown:
        print(f"Error: Unrecognised arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    map_id_out = os.path.join(args.output_dir, f"{args.output_prefix}.map_id.tsv")
    if os.path.exists(map_id_out):
        os.remove(map_id_out)

    # -- read reference maps --
    print("Reading VHH file ...", flush=True)
    sep  = detect_separator(args.vhh_id_file)
    cols = args.vhh_cols.split(",")
    vhh_ids = pl.read_csv(args.vhh_id_file, has_header=True, columns=cols, separator=sep)
    vhh_ids = vhh_ids.rename({cols[0]: "vhh_id", cols[1]: "vhh_seq"})
    vhh_map = dict(zip(vhh_ids["vhh_seq"], vhh_ids["vhh_id"]))
    del vhh_ids; gc.collect()

    print("Reading target file ...", flush=True)
    sep  = detect_separator(args.target_id_file)
    cols = args.target_cols.split(",")
    target_ids = pl.read_csv(args.target_id_file, has_header=True, columns=cols, separator=sep)
    target_ids = target_ids.rename({cols[0]: "target_id", cols[1]: "target_seq"})
    target_map = dict(zip(target_ids["target_seq"], target_ids["target_id"]))
    del target_ids; gc.collect()

    print("Mapping IDs to records ...", flush=True)
    list_results = []
    for i, chunk_result in enumerate(
        process_records_in_chunk(args.associate_file, args, vhh_map, target_map)
    ):
        print(
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} --> "
            f"Processed chunk {i+1} with {chunk_result.height} effective records",
            flush=True,
        )
        if not chunk_result.is_empty():
            list_results.append(chunk_result)

    print("Finished processing all records.", flush=True)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Writing output ...", flush=True)

    list_results_filtered = [df for df in list_results if df.height > 0]
    if list_results_filtered:
        df_records = pl.concat(list_results_filtered, how="vertical")
        df_records.write_csv(map_id_out, separator="\t")
    else:
        with open(map_id_out, "w") as f:
            f.write("something wrong with input files, no records found!\n")
        sys.exit(0)