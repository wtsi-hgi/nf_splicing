#-- import modules --#
import os
import sys
import argparse
import re
import csv
import polars as pl
import pyarrow
from datetime import datetime

#-- functions --#
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px != py:
            self.parent[py] = px

def parse_and_cluster(df: pl.DataFrame, tol: int) -> pl.DataFrame:
    """
    Parse junction bed file and cluster junctions based on donor/acceptor positions.
    Parameters:
        -- df: polars DataFrame, junction bed file
        -- tol: int, tolerance for clustering
    Returns:
        -- df: polars DataFrame, junction bed file with cluster ids
    """
    strands, donors, acceptors = [], [], []
    for strand, start, sizes, starts in zip(df["strand"], df["start"], df["blockSizes"], df["blockStarts"]):
        sizes = list(map(int, sizes.split(",")))
        starts = list(map(int, starts.split(",")))
        strands.append(strand)
        donors.append(start + sizes[0])
        acceptors.append(start + starts[1])

    df = df.with_columns([pl.Series("donor", donors), pl.Series("acceptor", acceptors)])

    n = df.height
    uf = UnionFind(n)

    # Compare all pairs of junctions to find clusters
    for i in range(n):
        for j in range(i + 1, n):
            if strands[i] != strands[j]:
                continue
            if abs(donors[i] - donors[j]) <= tol and abs(acceptors[i] - acceptors[j]) <= tol:
                uf.union(i, j)

    # Assign cluster ids based on union-find results
    clusters_ids = [uf.find(i) for i in range(n)]
    return df.with_columns(pl.Series("cluster_id", clusters_ids))

def merge_clusters(df: pl.DataFrame) -> pl.DataFrame:
    """
    Merge clustered junctions and summarize coverage.
    Parameters:
        -- df: polars DataFrame, junction bed file with cluster ids
    Returns:
        -- df: polars DataFrame, merged junctions with summarized coverage
    """
    merged = []
    n_groups = df.select(pl.struct(["chrom", "strand", "cluster_id"])).unique().height
    for i, ((chrom, strand, cluster_id), df_sub) in enumerate(df.group_by(["chrom", "strand", "cluster_id"]), start = 1):
        dominant_junction = df_sub.sort("coverage", descending = True).head(1)
        total_score = df_sub["coverage"].sum()

        merged.append(pl.DataFrame({ "chrom":    [chrom],
                                     "start":    dominant_junction["start"],
                                     "end":      dominant_junction["end"],
                                     "name":     dominant_junction['name'],
                                     "coverage": [total_score],
                                     "strand":   [strand],
                                     "donor":    dominant_junction["donor"],
                                     "acceptor": dominant_junction["acceptor"] }))
        
        if i % 10000 == 0:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed {i} / {n_groups} clusters", flush = True)

    df_merged = pl.concat(merged, how = "vertical")
    return df_merged.sort(["chrom", "donor", "acceptor"])

def annotate_junction(var_id, donor_pos, acceptor_pos, exon_pos, intron_pos, min_overlap):
    """
    Annotate junctions based on exon and intron positions.
    Parameters:
        -- var_id: str, variant id (chromosome)
        -- donor_pos: int, junction start position
        -- acceptor_pos: int, junction end position
        -- exon_pos: polars DataFrame, exon positions
        -- intron_pos: polars DataFrame, intron positions
        -- min_overlap: int, minimum overlap to consider partial splicing
    Returns:
        -- annotation: str, annotation string
    """
    var_id_base = var_id if args.lib_type in ["random_intron", "random_exon"] else var_id.split("/")[0]
    
    exon_pos_var = exon_pos.filter(pl.col("var_id") == var_id_base).to_pandas()
    intron_pos_var = intron_pos.filter(pl.col("var_id") == var_id_base).to_pandas()

    if exon_pos_var.empty:
        return "exon_data_missing"
    
    if intron_pos_var.empty:
        return "intron_data_missing"

    if donor_pos < exon_pos_var["exon_start"].min() or acceptor_pos > exon_pos_var["exon_end"].max():
        return "out_of_range"

    # intron retention
    for i, intron in intron_pos_var.iterrows():
        # intron retention: allow 2bp mismatch for start and end
        # because of 0-based start and 1-based end in bed format, how to define donor and acceptor positions
        check_start = (intron["intron_start"] - 2 <= donor_pos <= intron["intron_start"] + 2)
        check_end   = (intron["intron_end"] - 2   <= acceptor_pos   <= intron["intron_end"] + 2)
        if check_start and check_end:
            if i == 0:
                return f"intron_retention_{intron_pos_var.iloc[i+1]['intron_id']}"
            elif i > 0:
                return f"intron_retention_{intron_pos_var.iloc[i-1]['intron_id']}"

    # exon splicing
    annostr = []
    for _, exon in exon_pos_var.iterrows():
        if donor_pos > exon["exon_start"] + min_overlap and donor_pos < exon["exon_end"] - min_overlap:
            annostr.append(f"exon_splicing_3p_{exon['exon_id']}")

        if acceptor_pos > exon["exon_start"] + min_overlap and acceptor_pos < exon["exon_end"] - min_overlap:
            annostr.append(f"exon_splicing_5p_{exon['exon_id']}")

        if donor_pos < exon["exon_start"] + min_overlap and acceptor_pos > exon["exon_end"] - min_overlap:
            annostr.append(f"exon_skipping_{exon['exon_id']}")

    # intron splicing
    for idx, intron in intron_pos_var.iterrows():
        if donor_pos > intron["intron_start"] + min_overlap and donor_pos < intron["intron_end"] - min_overlap:
            annostr.append(f"intron_retension_5p_{intron['intron_id']}")
        if acceptor_pos > intron["intron_start"] + min_overlap and acceptor_pos < intron["intron_end"] - min_overlap:
            annostr.append(f"intron_retension_3p_{intron['intron_id']}")

    return ";".join(annostr) if annostr else "no_annotation"

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Classifying junctions types from a junction bed file.", allow_abbrev = False)
    parser.add_argument("--lib_type",            type = str, required = True,       help = "Library type (random_intron, random_exon, muta_intron, muta_exon)")
    parser.add_argument("--bed_file",            type = str, required = True,       help = "bed file generated from regtools junctions extract")
    parser.add_argument("--exon_pos",            type = str, required = True,       help = "exon position file")
    parser.add_argument("--cluster_tol",         type = int, default = 2,           help = "maximum tolerance of donor/acceptor positions for clustering")
    parser.add_argument("--junc_cov",            type = int, default = 2,           help = "minimum junction coverage to keep")
    parser.add_argument("--min_overlap",         type = int, default = 2,           help = "minimum anchor to consider partial splicing")
    parser.add_argument("--output_dir",          type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",       type = str, default = '',          help = "output prefix")

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Error: Unrecognized arguments: {' '.join(unknown)}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output_prefix == '':
        output_prefix = os.path.splitext(os.path.basename(args.bed_file))[0]
    else:
        output_prefix = args.output_prefix

    # -- read input files -- #
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Reading bed file, please wait...", flush = True)
    df_junctions = pl.read_csv(args.bed_file, separator = "\t", has_header = False)
    df_junctions = df_junctions.rename({ "column_1": "chrom",
                                         "column_2": "start",
                                         "column_3": "end",
                                         "column_4": "name",
                                         "column_5": "coverage",
                                         "column_6": "strand",
                                         "column_11": "blockSizes",
                                         "column_12": "blockStarts" })

    df_exon_pos = pl.read_csv(args.exon_pos, separator = "\t", has_header = False)
    df_exon_pos = df_exon_pos.rename({ "column_1": "var_id",
                                       "column_2": "exon_id",
                                       "column_3": "exon_start",
                                       "column_4": "exon_end" })

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    junction_out = f"{output_prefix}.classified_junctions.tsv"
    if os.path.exists(junction_out):
        os.remove(junction_out)
    
    # -- processing -- #
    chroms = df_junctions["chrom"].unique().to_list()
    total_chroms = len(chroms)
    grouped_clusters = []
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing bed file, please wait...", flush = True)
    for chrom, df_sub in df_junctions.group_by("chrom"):
        df_sub_clusters = parse_and_cluster(df_sub, args.cluster_tol)
        grouped_clusters.append(df_sub_clusters)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Total number of variants: {total_chroms}", flush = True)
    
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Merging splicing junctions by clusters, please wait...", flush = True)
    df_junctions_clusters = pl.concat(grouped_clusters, how = "vertical")
    df_junctions_merged = merge_clusters(df_junctions_clusters)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Total number of splicing junctions: {df_junctions_merged.height}", flush = True)

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Classify splicing junctions, please wait...", flush = True)
    df_intron_pos = ( df_exon_pos.with_columns(prev_end = pl.col("exon_end").shift().over("var_id"),
                                               prev_chrom = pl.col("var_id").shift().over("var_id"))
                                 .filter(pl.col("var_id") == pl.col("prev_chrom"))
                                 .with_columns([(pl.col("prev_end") + 1).alias("intron_start"),
                                                (pl.col("exon_start") - 1).alias("intron_end")])
                                 .select(["var_id", "intron_start", "intron_end"])
                                 .with_columns([pl.format("I{}", pl.arange(1, pl.len() + 1).over("var_id")).alias("intron_id")]))

    df_junctions_filtered = df_junctions_merged.filter(pl.col("coverage") >= args.junc_cov)
    df_junctions_classes = df_junctions_filtered.with_columns([
        pl.struct(["chrom", "donor", "acceptor"]).map_elements(
            lambda row: annotate_junction(row["chrom"], row["donor"], row["acceptor"], df_exon_pos, df_intron_pos, args.min_overlap),
            return_dtype=pl.String
        ).alias("annotation")
    ])
    df_junctions_classes = df_junctions_classes.filter(pl.col("annotation") != "no_annotation")
    df_junctions_classes.write_csv(junction_out, separator = "\t")

    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Finished", flush = True)