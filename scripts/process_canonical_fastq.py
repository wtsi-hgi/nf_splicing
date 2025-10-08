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
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
from sequence_utils import pigz_open, fastq_iter, read_first_fasta_seq, match_approximate, extract_sequence, reverse_complement, check_barcode

def process_se_read(read: tuple) -> list:
    """
    Detect canonical splicing event in a read
    Parameters:
        -- read: tuple (header, sequence, quality) for read
    Returns:
        -- tuple: (variant_seq, barcode_seq) extracted sequences
    """
    read_seq = read[1]
    dict_barcode = { 'read_ref'    : "NA",
                     'var_id'      : "NA",
                     'barcode'     : "NA",
                     'ref_type'    : "NA" }

    barcode_seq = extract_sequence(read_seq, args.barcode_up, args.barcode_down, args.max_mismatch)
    if (barcode_seq in {"upstream not found", "downstream not found"}):
        barcode_seq = extract_sequence(reverse_complement(read_seq), args.barcode_up, args.barcode_down, args.max_mismatch)
        if (barcode_seq in {"upstream not found", "downstream not found"}):
            return read, {}
        else:
            read_seq = reverse_complement(read[1])

    if args.barcode_check:
        check_res = check_barcode(barcode_seq, args.barcode_temp, args.barcode_mismatch)
        if check_res is None:
            return read, {}
        else:
            barcode_seq = check_res

    dict_barcode['barcode'] = barcode_seq

    res_bar_var = dict_bar_var.get(barcode_seq)
    if res_bar_var is not None:
        dict_barcode['read_ref'] = res_bar_var["var_id"]
        dict_barcode['var_id'] = res_bar_var["var_id"]
        exon = res_bar_var["exon"]
        ref_canonical_inclusion = pre_exon + exon + post_exon
        ref_canonical_skipping = pre_exon + post_exon
        ref_found = match_approximate(read_seq, ref_canonical_inclusion, args.max_mismatch)
        if ref_found == -1:
            ref_found = match_approximate(read_seq, ref_canonical_skipping, args.max_mismatch)
            if ref_found == -1:
                return read, {}
            else:
                dict_barcode['ref_type'] = "canonical_skipping"
                return (), dict_barcode
        else:
            dict_barcode['ref_type'] = "canonical_inclusion"
            return (), dict_barcode
    else:
        return read, {}

def batch_process_se_reads(batch_reads: list) -> list:
    """
    Process a batch of read pairs to extract variants and barcodes.
    Parameters:
        -- batch_reads: list of tuples, each containing (read1, read2)
    Returns:
        -- list of tuples: each tuple contains (variant_seq, barcode_seq) for processed read
    """
    results = []
    for read in batch_reads:
        result = process_se_read(read)
        results.append(result)
    return results

def function_processpool_se(args):
    """
    Wrapper function for process pool as ProcessPoolExecutor expects a function rather than returned results.
    """
    return batch_process_se_reads(args)

def process_se_reads_in_chunk(read_path):
    """
    Read paired-end FASTQ files in chunks and process each pair of reads in parallel
    Parameters:
        -- read_path: path to single-end read FASTQ file
    Yields:
        -- list of tuples: each tuple contains (variant_seq, barcode_seq) for processed read
    """
    read_fh = io.TextIOWrapper(pigz_open(read_path).stdout) if read_path.endswith(".gz") else open(read_path)
    read_iter = fastq_iter(read_fh)

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        while True:
            read_chunk = list(islice(read_iter, args.chunk_size))
            if not read_chunk:
                break

            # Divide chunk into batches
            # if process_long_read is very fast, we can use a larger batch size to make better use of CPU resources
            # if process_long_read is very slow, we can use a smaller batch size to make better use of CPU resources
            # batch_size = min(args.chunk_size, 5000)
            batch_size = min(args.chunk_size, 100)
            read_batches = [
                read_chunk[i:i+batch_size]
                for i in range(0, len(read_chunk), batch_size)
            ]

            args_list = [ batch for batch in read_batches ]
            batch_results = list(executor.map(function_processpool_se, args_list))

            flat_results = [item for batch in batch_results for item in batch]
            fail_reads, pass_barcodes = zip(*flat_results)

            for read in fail_reads:
                if read:
                    fastq_fh.write(f"{read[0]}\n{read[1]}\n+\n{read[2]}\n".encode("utf-8"))

            filtered_pass_barcodes = [dict for dict in pass_barcodes if dict]
            yield pl.DataFrame(filtered_pass_barcodes)

    read_fh.close()

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
        output_prefix = os.path.splitext(os.path.basename(args.bam_file))[0]
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
                     for row in df_bar_var.iter_rows(named=True) }

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
    else:
        fail_r1_fastq = f"{output_prefix}.filter_se.fail.r1.fastq.gz"
        fail_r2_fastq = f"{output_prefix}.filter_se.fail.r2.fastq.gz"

    barcode_out = f"{output_prefix}.canonical_barcodes.tsv"
    if os.path.exists(barcode_out):
        os.remove(barcode_out)
    
    # -- parallel processing -- #
    barcode_list = []
    pigz_proc = subprocess.Popen( ["pigz", "-c", "-p", str(args.threads)],
                                  stdin=subprocess.PIPE,
                                  stdout=open(fail_fastq, "wb") )
    if args.read_type == 'se':
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing fastq file, please wait...", flush = True)
        with pigz_proc.stdin as fastq_fh:
            for i, chunk_result in enumerate(process_se_reads_in_chunk(path_read)):
                print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed chunk {i+1} with {args.chunk_size} reads", flush = True)
                barcode_list.append(chunk_result)
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Finished.", flush = True)

        pigz_proc.wait()

    # -- clean and format the extracted barcodes from reads -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Generating barcode results, please wait...", flush = True)
    filtered_barcode_list = [df for df in barcode_list if df.shape[1] > 0]
    if filtered_barcode_list:
        df_barcode = pl.concat(filtered_barcode_list, how = "vertical")
        df_barcode_counts = df_barcode.group_by(["read_ref", "var_id", "barcode", "ref_type"]).agg(pl.len().alias("count"))
        df_barcode_counts = df_barcode_counts.select(["read_ref", "var_id", "barcode", "count", "ref_type"])
        df_barcode_counts.write_csv(barcode_out, separator = "\t", null_value = "NA")
    else:
        with open(barcode_out, "w") as f:
            f.write("no barcode found in the reads, please check your barcode marker or template!\n")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Finished.", flush = True)
