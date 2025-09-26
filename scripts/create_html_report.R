#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "vroom", "ggVennDiagram", "htmltools", "reactable", "optparse", "sparkline", "UpSetR", "patchwork", "glue", "scales", "ggExtra")
invisible(lapply(packages, quiet_library))

# -- options -- #
option_list <- list(make_option(c("-r", "--rscript_dir"),          type = "character", help = "directory path of R scripts",               default = NULL),
                    make_option(c("-l", "--lib_type"),             type = "character", help = "library type",                              default = NULL),
                    make_option(c("-e", "--exon_pos"),             type = "character", help = "exon position file",                        default = NULL),
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
if(is.null(opt$lib_type))             stop("-l, library type is required!", call. = FALSE)
if(is.null(opt$barcode_association))  stop("-b, barcode association file is required!", call. = FALSE)
if(is.null(opt$exon_pos))             stop("-e, exon position file is required!", call. = FALSE)
if(is.null(opt$sample_id))            stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$trim_stats))           stop("-t, list of trim stats files is required!", call. = FALSE)
if(is.null(opt$merge_stats))          stop("-m, list of merge stats files is required!", call. = FALSE)
if(is.null(opt$bwa_idxstats))         stop("-f, list of bwa map idxstats files is required!", call. = FALSE)
if(is.null(opt$hisat2_stats))         stop("-a, list of hisat2 map summary files is required!", call. = FALSE)
if(is.null(opt$canonical_barcodes))   stop("-c, list of extracted canonical barcode files is required!", call. = FALSE)
if(is.null(opt$novel_barcodes))       stop("-n, list of extracted novel barcode files is required!", call. = FALSE)
if(is.null(opt$classified_junctions)) stop("-j, list of classified junction file is required!", call. = FALSE)
if(is.null(opt$splicing_counts))      stop("-d, list of splicing counts is required!", call. = FALSE)

valid_lib_types <- c("random_intron", "random_exon", "muta_intron", "muta_exon")
if(!(opt$lib_type %in% valid_lib_types))
{
    stop(paste0("-l, library type must be one of: ", paste(valid_lib_types, collapse = ", ")), call. = FALSE)
}

# -- modules -- #
source(file.path(opt$rscript_dir, "report_utils.R"))
source(file.path(opt$rscript_dir, "report_plots.R"))
source(file.path(opt$rscript_dir, "report_html.R"))

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
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Reading input files ...")

barcode_association <- as.data.table(vroom(opt$barcode_association, delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
exon_pos <- as.data.table(vroom(opt$exon_pos, delim = "\t", comment = "#", col_names = FALSE, show_col_types = FALSE))
setnames(exon_pos, c("var_id", "exon_id", "exon_start", "exon_end"))

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
    tmp_value <- as.numeric(str_extract(grep("reads passed filter:", readLines(files_trim_stats[i]), value = TRUE), "\\d+"))
    total_reads[i] <- tmp_value / 2

    merged_reads[i] <- as.numeric(str_extract(grep("Combined pairs:", readLines(files_merge_stats[i]), value = TRUE), "\\d+"))
    unmerged_reads[i] <- as.numeric(str_extract(grep("Uncombined pairs:", readLines(files_merge_stats[i]), value = TRUE), "\\d+"))

    inclusion_reads[i] <- as.numeric(read.table(files_bwa_idxstats[i])[1, 2])
    skipping_reads[i] <- as.numeric(read.table(files_bwa_idxstats[i])[2, 2])

    map_reads[i] <- as.numeric(read.table(files_hisat2_stats[i])[1, 2])
    unexplain_reads[i] <- as.numeric(read.table(files_hisat2_stats[i])[2, 2]) - as.numeric(read.table(files_hisat2_stats[i])[1, 2])

    canonical_barcodes[[sample_reps[i]]] <- as.data.table(vroom(files_canonical_barcodes[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    novel_barcodes[[sample_reps[i]]] <- as.data.table(vroom(files_novel_barcodes[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

    classified_junctions[[sample_reps[i]]] <- as.data.table(vroom(files_classified_junctions[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    splicing_counts[[sample_reps[i]]] <- as.data.table(vroom(files_splicing_counts[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
}

# -- processing -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Preparing tables and figures ...")

# 1. bar plots for statistics of sequencing reads
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> Creating read summaries and plots ...")

data_barplots <- create_barplots(barcode_association, sample_reps, total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads)
summary_reads <- data_barplots[[1]]
summary_pct <- data_barplots[[2]]
barplots_combined <- wrap_plots(data_barplots[[3]], nrow = 1)

fwrite(summary_reads, file = paste0(sample_prefix, ".summary_reads.txt"), sep = "\t", row.names = FALSE)
fwrite(summary_pct, file = paste0(sample_prefix, ".summary_pct.txt"), sep = "\t", row.names = FALSE)

png(paste0(sample_prefix, ".reads_pct.png"), width = 800, height = 800, units = "px", res = 130)
barplots_combined
dev.off()

# 2. venn diagrams for detected barcodes and talbes of detected barcodes and variants
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> Creating barcode summaries and plots ...")

venn_diagrams <- create_venn_diagrams(canonical_barcodes, novel_barcodes, barcode_association)
venn_diagrams_combined <- wrap_plots(venn_diagrams, nrow = 1)

png(paste0(sample_prefix, ".barcodes_venn.png"), width = 1800, height = 600, units = "px", res = 120)
venn_diagrams_combined
dev.off()

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

# 3. venn diagrams for detected junctions and correlation of junction quantifications
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> Creating junction correlation plot ...")

junction_plots <- create_junction_plots(classified_junctions)

png(paste0(sample_prefix, ".junctions_venn.png"), width = 800, height = 800, units = "px", res = 150)
print(junction_plots[[1]])
dev.off()

junctions_category <- junction_plots[[2]]
fwrite(junctions_category, file = paste0(sample_prefix, ".junctions_category.txt"), sep = "\t", row.names = FALSE)

cor_data <- junctions_category[, ..sample_reps]
cor_data[is.na(cor_data)] <- 0
cor_data_norm <- cor_data %>% mutate(across(everything(), ~ log2((. / sum(.)) * 1e6 + 1)))

png(paste0(sample_prefix, ".junctions_corr.png"), width = 1200, height = 1200, units = "px", res = 150)
pairs(cor_data_norm,
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = function(x, y, ...) {panel.smooth(x, y, method = "lm", ...)},
      use = "complete.obs")
dev.off()

# 4. upset plots for splicing events and tables of splicing events
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> Creating junction category plot ...")

junctions_category_reshape <- junctions_category %>%
                                separate_rows(annotation, sep = ";") %>%
                                select(c(var_id, annotation, avg_cov)) %>%
                                group_by(var_id, annotation) %>%
                                summarise(avg_cov = sum(avg_cov, na.rm = TRUE), .groups = "drop") %>%
                                mutate(log_avg_cov = log2(avg_cov + 1)) %>%
                                select(c(var_id, annotation, log_avg_cov)) %>% 
                                filter(annotation != "no_annotation")
upset_input <- junctions_category_reshape %>%
                distinct(var_id, annotation) %>%
                mutate(value = 1) %>%
                pivot_wider(names_from = annotation, values_from = value, values_fill = 0)
upset_join <- junctions_category_reshape %>% left_join(upset_input, by = "var_id")

png(paste0(sample_prefix, ".junctions_category.png"), width = 2400, height = 1600, units = "px", res = 200)
upset(as.data.frame(upset_join), 
      sets = unique(upset_join$annotation), 
      order.by = "freq", 
      matrix.color = "yellowgreen",
      main.bar.color = "royalblue", 
      sets.bar.color = "yellowgreen",
      boxplot.summary = c("log_avg_cov"))
dev.off()

# 5. psi correlation of splicing events
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> Creating PSI plot ...")

for(i in 1:length(sample_reps))
{
    splicing_counts[[sample_reps[i]]] <- splicing_counts[[sample_reps[i]]] %>%
                                            rowwise() %>%
                                            mutate(PSI = canonical_inclusion / sum(c_across(-var_id))) %>%
                                            ungroup()
}

psi_list <- lapply(names(splicing_counts), function(rep_id) { splicing_counts[[rep_id]] %>%
                                                                select(var_id, PSI) %>%
                                                                rename(!!rep_id := PSI) })
psi_wide <- reduce(psi_list, full_join, by = "var_id")
fwrite(psi_wide, file = paste0(sample_prefix, ".psi_values.txt"), sep = "\t", row.names = FALSE)

png(paste0(sample_prefix, ".psi_correlation.png"), width = 1200, height = 1200, units = "px", res = 100)
pairs(psi_wide[, -1],
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = function(x, y, ...) {panel.smooth(x, y, method = "lm", ...)},
      use = "complete.obs")
dev.off()

# 6. junction distribution plots
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |--> Creating junction distribution plot ...")

data_rescale <- rescale_junctions(junctions_category, exon_pos, opt$lib_type)
junctions_rescale <- data_rescale[[1]]
exons_template <- data_rescale[[2]]
introns_template <- data_rescale[[3]]

junctions_distri <- create_junction_distribution(junctions_rescale, exons_template, introns_template)
junctions_range <- junctions_distri[[1]]
junctions_diagramplot <- junctions_distri[[2]]
junctions_scatterplot <- junctions_distri[[3]]

for(i in seq_along(junctions_range))
{
    png(paste0(sample_prefix, ".junctions_diagram_range_", i, ".png"), width = 1600, height = 600, units = "px", res = 200)
    print(junctions_diagramplot[[i]])
    dev.off()

    png(paste0(sample_prefix, ".junctions_scatter_range_", i, ".png"), width = 1600, height = 1400, units = "px", res = 200)
    print(junctions_scatterplot[[i]])
    dev.off()
}

# -- reporting -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Creating final html report...")

file_summary_reads <- paste0(sample_prefix, ".summary_reads.txt")
file_summary_pct <- paste0(sample_prefix, ".summary_pct.txt")
plot_reads_pct <- paste0(sample_prefix, ".reads_pct.png")
plot_barcodes_venn <- paste0(sample_prefix, ".barcodes_venn.png")
file_barcodes_summary <- paste0(sample_prefix, ".summary_barvars.txt")
plot_junctions_venn <- paste0(sample_prefix, ".junctions_venn.png")
plot_junctions_corr <- paste0(sample_prefix, ".junctions_corr.png")
plot_junctions_category <- paste0(sample_prefix, ".junctions_category.png")
file_junctions_category <- paste0(sample_prefix, ".junctions_category.txt")
list_files_junctions_diagram <- list.files(pattern = paste0(sample_prefix, ".junctions_diagram_range_.*.png$"))
names(list_files_junctions_diagram) <- names(junctions_range)
list_files_junctions_scatter <- list.files(pattern = paste0(sample_prefix, ".junctions_scatter_range_.*.png$"))
names(list_files_junctions_scatter) <- names(junctions_range)
plot_psi_corr <- paste0(sample_prefix, ".psi_correlation.png")
file_psi <- paste0(sample_prefix, ".psi_values.txt")

file_render_context <- paste0(sample_prefix, ".html_report.Rmd")
create_html_render(file_summary_reads, 
                   file_summary_pct,
                   plot_reads_pct,
                   plot_barcodes_venn,
                   file_barcodes_summary,
                   plot_junctions_venn,
                   plot_junctions_corr,
                   plot_junctions_category,
                   file_junctions_category
                   list_files_junctions_diagram,
                   list_files_junctions_scatter,
                   plot_psi_corr,
                   file_psi,
                   file_render_context)

rmarkdown::render(file_render_context, clean = TRUE, quiet = TRUE)
invisible(file.remove(file_render_context))
