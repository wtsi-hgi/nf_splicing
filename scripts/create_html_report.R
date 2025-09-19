#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "vroom", "ggVennDiagram", "htmltools", "reactable", "optparse", "sparkline", "UpSetR", "patchwork", "glue")
invisible(lapply(packages, quiet_library))

# -- options -- #
option_list <- list(make_option(c("-r", "--rscript_dir"),          type = "character", help = "directory path of R scripts",               default = NULL),
                    make_option(c("-b", "--barcode_association"),  type = "character", help = "barcode association file",                  default = NULL),
                    make_option(c("-s", "--sample_id"),            type = "character", help = "list of sample IDs",                        default = NULL),
                    make_option(c("-t", "--trim_stats"),           type = "character", help = "list of trim stats files",                  default = NULL),
                    make_option(c("-m", "--merge_stats"),          type = "character", help = "list of merge stats files",                 default = NULL),
                    make_option(c("-f", "--bwa_idxstats"),         type = "character", help = "list of bwa map idxstats files",            default = NULL),
                    make_option(c("-a", "--hisat2_stats"),         type = "character", help = "list of hisat2 map summary files",          default = NULL),
                    make_option(c("-c", "--canonical_barcodes"),   type = "character", help = "list of extracted canonical barcode files", default = NULL),
                    make_option(c("-n", "--novel_barcodes"),       type = "character", help = "list of extracted novel barcode files",     default = NULL),
                    make_option(c("-j", "--classified_junctions"), type = "character", help = "list of classified junction files",         default = NULL),
                    make_option(c("-d", "--splicing_counts"),      type = "character", help = "list of splicing counts",                   default = NULL),
                    make_option(c("-o", "--output_dir"),           type = "character", help = "output directory",                          default = getwd()),
                    make_option(c("-p", "--prefix"),               type = "character", help = "output prefix",                             default = "sample"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

# -- check options -- #
if(is.null(opt$rscript_dir))          stop("-r, directory path of R scripts is required!", call. = FALSE)
if(is.null(opt$barcode_association))  stop("-b, barcode association file is required!", call. = FALSE)
if(is.null(opt$sample_id))            stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$trim_stats))           stop("-t, list of trim stats files is required!", call. = FALSE)
if(is.null(opt$merge_stats))          stop("-m, list of merge stats files is required!", call. = FALSE)
if(is.null(opt$bwa_idxstats))         stop("-f, list of bwa map idxstats files is required!", call. = FALSE)
if(is.null(opt$hisat2_stats))         stop("-a, list of hisat2 map summary files is required!", call. = FALSE)
if(is.null(opt$canonical_barcodes))   stop("-c, list of extracted canonical barcode files is required!", call. = FALSE)
if(is.null(opt$novel_barcodes))       stop("-n, list of extracted novel barcode files is required!", call. = FALSE)
if(is.null(opt$classified_junctions)) stop("-j, list of classified junction file is required!", call. = FALSE)
if(is.null(opt$splicing_counts))      stop("-d, list of splicing counts is required!", call. = FALSE)

# -- modules -- #
source(file.path(opt$rscript_dir, "report_utils.R"))
source(file.path(opt$rscript_dir, "report_utils.R"))

# -- inputs -- #
sample_reps                <- unlist(strsplit(opt$sample_id, ","))
files_trim_stats           <- unlist(strsplit(opt$trim_stats, ","))
files_merge_stats          <- unlist(strsplit(opt$merge_stats, ","))
files_bwa_idxstats         <- unlist(strsplit(opt$bwa_idxstats, ","))
files_hisat2_stats         <- unlist(strsplit(opt$hisat2_stats, ","))
files_canonical_barcodes   <- unlist(strsplit(opt$canonical_barcodes, ","))
files_novel_barcodes       <- unlist(strsplit(opt$novel_barcodes, ","))
files_classified_junctions <- unlist(strsplit(opt$classified_junctions, ","))
files_splicing_counts      <- unlist(strsplit(opt$splicing_counts, ","))

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

output_prefix <- opt$prefix

# -- reading files -- #
barcode_association <- as.data.table(vroom(opt$barcode_association, delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

total_reads <- numeric(length(sample_reps))

merged_reads <- numeric(length(sample_reps))
unmerged_reads <- numeric(length(sample_reps))

inclusion_reads <- numeric(length(sample_reps))
skipping_reads <- numeric(length(sample_reps))

map_reads <- numeric(length(sample_reps))
unexplain_reads <- numeric(length(sample_reps))

canonical_barcodes <- list()
novel_barcodes <- list()

classified_junctions <- list()
splicing_counts <- list()

for(i in seq_along(sample_reps))
{
    tmp_value <- as.numeric(str_extract(grep("reads passed filter:", readLines(trim_files[i]), value = TRUE), "\\d+"))
    total_reads[i] <- tmp_value / 2

    merged_reads[i] <- as.numeric(str_extract(grep("Combined pairs:", readLines(merge_files[i]), value = TRUE), "\\d+"))
    unmerged_reads[i] <- as.numeric(str_extract(grep("Uncombined pairs:", readLines(merge_files[i]), value = TRUE), "\\d+"))

    inclusion_reads[i] <- as.numeric(read.table(filter_files[i])[1, 2])
    skipping_reads[i] <- as.numeric(read.table(filter_files[i])[2, 2])

    map_reads[i] <- as.numeric(read.table(map_files[i])[1, 2])
    unexplain_reads[i] <- as.numeric(read.table(map_files[i])[2, 2]) - as.numeric(read.table(map_files[i])[1, 2])

    canonical_barcodes[[reps[i]]] <- as.data.table(vroom(canonical_barcode_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    novel_barcodes[[reps[i]]] <- as.data.table(vroom(novel_barcode_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

    classified_junctions[[reps[i]]] <- as.data.table(vroom(classified_junctions_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    splicing_counts[[reps[i]]] <- as.data.table(vroom(splicing_counts_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
}

# -- processing -- #
barplots <- create_barplots(sample_reps, total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads)
barplots_combined <- wrap_plots(barplots, nrow = 1)

venn_diagrams <- create_venn_diagrams(canonical_barcodes, novel_barcodes, barcode_association)
venn_diagrams_combined <- wrap_plots(venn_diagrams, nrow = 1)

num_detected_barcodes <- numeric(length(sample_reps))
num_detected_variants <- numeric(length(sample_reps))
for(i in seq_along(sample_reps))
{
    sequencing_barcodes <- unique(c(canonical_barcodes[[i]]$barcode, novel_barcodes[[i]]$barcode))
    num_detected_barcodes[i] <- length(intersect(sequencing_barcodes, barcode_association$barcode))

    sequencing_variants <- unique(c(canonical_barcodes[[i]]$var_id, novel_barcodes[[i]]$var_id))
    num_detected_variants[i] <- length(sequencing_variants) - 1 # there is NAs in var_id
}

summary_barvars <- as.data.table(cbind(rep(length(barcode_association$barcode), 3),
                                       rep(length(unique(barcode_association$var_id)), 3),
                                       num_detected_barcodes,
                                       num_detected_variants))
colnames(summary_barvars) <- c("num_expected_barcodes", "num_expected_variants", "num_detected_barcodes", "num_detected_variants")
summary_barvars <- summary_barvars %>% 
                        mutate(pct_detected_barcodes = 100 * num_detected_barcodes / num_expected_barcodes,
                               pct_detected_variants = 100 * num_detected_variants / num_expected_variants)



