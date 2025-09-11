#-- import modules --#
import os
import sys
import argparse
import re
import csv
import polars as pl
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
    for (chrom, strand, cluster_id), df_sub in df.group_by(["chrom", "strand", "cluster_id"]):
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

    df_merged = pl.concat(merged, how = "vertical")
    return df_merged.sort(["chrom", "donor", "acceptor"])

#-- main execution --#
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Classifying junctions types from a junction bed file.", allow_abbrev = False)
    parser.add_argument("--bed_file",            type = str, required = True,       help = "bed file generated from regtools junctions extract")
    parser.add_argument("--exon_pos",            type = str, required = True,       help = "exon position file")
    parser.add_argument("--cluster_tol",         type = int, default = 2,           help = "maximum tolerance of donor/acceptor positions for clustering")
    parser.add_argument("--junc_cov",            type = int, default = 2,           help = "minimum junction coverage to keep")
    parser.add_argument("--output_dir",          type = str, default = os.getcwd(), help = "output directory")
    parser.add_argument("--output_prefix",       type = str, default = '',          help = "output prefix")

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

    # -- prepare output files -- #
    os.makedirs(args.output_dir, exist_ok = True)
    os.chdir(args.output_dir)

    junction_out = f"{output_prefix}.junction_clusters.txt"
    if os.path.exists(junction_out):
        os.remove(junction_out)
    
    # -- processing -- #
    chroms = df_junctions["chrom"].unique().to_list()
    total_chroms = len(chroms)
    grouped_clusters = []
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Processing bed file, please wait...", flush = True)
    for i, (chrom, df_sub) in enumerate(df_junctions.group_by("chrom"), start = 1):
        df_sub_clusters = parse_and_cluster(df_sub, args.cluster_tol)
        grouped_clusters.append(df_sub_clusters)
        if i % 1000 == 0 or i == total_chroms:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} |--> Processed {i}/{total_chroms} variants...", flush = True)
    
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Merging splicing junctions by clusters, please wait...", flush = True)
    df_junctions_clusters = pl.concat(grouped_clusters, how = "vertical")
    df_junctions_merged = merge_clusters(df_junctions_clusters)
    # df_junctions_merged.write_csv(junction_out, separator = "\t")

    df_junctions_filtered = df_junctions_merged.filter(pl.col("coverage") >= args.junc_cov)

