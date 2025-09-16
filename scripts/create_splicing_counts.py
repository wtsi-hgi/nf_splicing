#-- import modules --#
import os
import sys
import argparse
import re
import csv
import polars as pl
from datetime import datetime


#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Creating the count matrix of splicing events.", allow_abbrev = False)
    parser.add_argument("--canonical_file", type = str, required = True,       help = "the input file of canonical counts")
    parser.add_argument("--novel_file",     type = str, required = True,       help = "the input file of novel counts")
    parser.add_argument("--barcode_file",   type = str, required = True,       help = "barcode association file")
    parser.add_argument("--output_dir",     type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",  type = str, default = '',          help = "output prefix")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output_prefix == '':
        output_prefix = os.path.splitext(os.path.basename(args.bam_file))[0]
    else:
        output_prefix = args.output_prefix
    
    # -- read input files -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Reading input files, please wait...", flush = True)
    df_canonical = pl.read_csv(args.canonical_file, separator = "\t", has_header = True)
    df_novel = pl.read_csv(args.novel_file, separator = "\t", has_header = True)
    df_bar_var = pl.read_csv(args.barcode_file, separator = "\t", has_header = True, columns = "var_id")

    categories = [ "intron_retention_I1",
                   "intron_retention_I2",
                   "intron_retension_5p_I1",
                   "intron_retension_5p_I2",
                   "intron_retension_3p_I1",
                   "intron_retension_3p_I2",
                   "exon_splicing_3p_E1",
                   "exon_splicing_3p_E2",
                   "exon_splicing_3p_E3",
                   "exon_splicing_5p_E1",
                   "exon_splicing_5p_E2",
                   "exon_splicing_5p_E3",
                   "exon_skipping_E2" ]
    
    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    counts_out = f"{output_prefix}.splicing_counts.txt"
    if os.path.exists(counts_out):
        os.remove(counts_out)

    # -- format data frame -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Formating data frame, please wait...", flush = True)
    df_canonical_format = (df_canonical.filter(pl.col("var_id") != "NA")
                                       .select(["var_id", "count", "ref_type"])
                                       .group_by(["var_id", "ref_type"])
                                       .agg(pl.col("count").sum()))
    df_canonical_pivot = df_canonical_format.pivot(values = "count", index = "var_id", on = "ref_type").fill_null(0)

    df_novel_format = (df_novel.select(["chrom", "coverage", "annotation"])
                               .with_columns(pl.col("annotation").str.split(";").alias("annotation_split"))
                               .drop("annotation")
                               .explode("annotation_split")
                               .filter(pl.col("annotation_split").is_in(categories))
                               .rename({"chrom": "var_id", "annotation_split": "annotation"})
                               .group_by(["var_id", "annotation"])
                               .agg(pl.col("coverage").sum()))
    df_novel_pivot = df_novel_format.pivot(values = "coverage", index = "var_id", on = "annotation").fill_null(0)
    missing_categories = [cat for cat in categories if cat not in df_novel_pivot.columns]
    if missing_categories:
        df_novel_pivot = df_novel_pivot.with_columns([pl.lit(0).alias(cat) for cat in missing_categories])
    df_novel_pivot = df_novel_pivot.select(["var_id"] + categories)

    # -- integrate the counts -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Integrating counts, please wait...", flush = True)
    df_bar_var = df_bar_var.unique(subset=["var_id"])
    df_bar_var = (df_bar_var.with_columns(pl.col("var_id").str.extract(r"(\d+)").cast(pl.Int64).alias("var_num"))
                            .sort("var_num")
                            .drop("var_num"))
    
    df_integrated = (df_bar_var.join(df_canonical_pivot, on = "var_id", how = "left")
                               .join(df_novel_pivot, on = "var_id", how = "left")
                               .fill_null(0))

    df_integrated.write_csv(counts_out, separator = "\t", null_value = "NA")
