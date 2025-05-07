#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "vroom", "ggVennDiagram", "htmltools", "reactable", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--barcode_association"), type = "character", help = "barcode association file",                  default = NULL),
                    make_option(c("-s", "--sample_id"),           type = "character", help = "list of sample IDs",                        default = NULL),
                    make_option(c("-t", "--trim_stats"),          type = "character", help = "list of trim stats files",                  default = NULL),
                    make_option(c("-m", "--merge_stats"),         type = "character", help = "list of merge stats files",                 default = NULL),
                    make_option(c("-f", "--filter_idxstats"),     type = "character", help = "list of bwa map idxstats files",            default = NULL),
                    make_option(c("-a", "--map_stats"),           type = "character", help = "list of hisat2 map summary files",          default = NULL),
                    make_option(c("-c", "--canonical_barcodes"),  type = "character", help = "list of extracted canonical barcode files", default = NULL),
                    make_option(c("-n", "--novel_barcodes"),      type = "character", help = "list of extracted novel barcode files",     default = NULL),
                    make_option(c("-o", "--output_dir"),          type = "character", help = "output directory",                          default = getwd()),
                    make_option(c("-p", "--prefix"),              type = "character", help = "output prefix",                             default = "sample"))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
  print_help(opt_parser)
  quit(status = 1)
}

# Check if required arguments are provided
if(is.null(opt$barcode_association)) stop("-b, barcode association file is required!", call. = FALSE)
if(is.null(opt$sample_id))           stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$trim_stats))          stop("-t, list of trim stats files is required!", call. = FALSE)
if(is.null(opt$merge_stats))         stop("-m, list of merge stats files is required!", call. = FALSE)
if(is.null(opt$filter_idxstats))     stop("-f, list of bwa map idxstats files is required!", call. = FALSE)
if(is.null(opt$map_stats))           stop("-a, list of hisat2 map summary files is required!", call. = FALSE)
if(is.null(opt$canonical_barcodes))  stop("-c, list of extracted canonical barcode files is required!", call. = FALSE)
if(is.null(opt$novel_barcodes))      stop("-n, list of extracted novel barcode files is required!", call. = FALSE)

#-- function --#
t_col <- function(col, rate)
{
    newcol <- rgb(col2rgb(col)["red", ],
                  col2rgb(col)["green", ],
                  col2rgb(col)["blue", ],
                  as.integer(rate * 255),
                  maxColorValue = 255)
    return(newcol)
}

select_colorblind <- function(col_id)
{
    col8 <- c("#D55E00", "#56B4E9", "#E69F00",
              "#009E73", "#F0E442", "#0072B2",
              "#CC79A7", "#000000")

    col12 <- c("#88CCEE", "#CC6677", "#DDCC77",
               "#117733", "#332288", "#AA4499",
               "#44AA99", "#999933", "#882255",
               "#661100", "#6699CC", "#888888")

    col15 <- c("red",       "royalblue", "olivedrab",
               "purple",    "violet",    "maroon1",
               "seagreen1", "navy",      "pink",
               "coral",     "steelblue", "turquoise1",
               "red4",      "skyblue",   "yellowgreen")

    col21 <- c("#F60239", "#009503", "#FFDC3D",
               "#9900E6", "#009FFA", "#FF92FD",
               "#65019F", "#FF6E3A", "#005A01",
               "#00E5F8", "#DA00FD", "#AFFF2A",
               "#00F407", "#00489E", "#0079FA",
               "#560133", "#EF0096", "#000000",
               "#005745", "#00AF8E", "#00EBC1")

    if (col_id == "col8") {
        return(col8)
    } else if (col_id == "col12") {
        return(col12)
    } else if (col_id == "col15") {
        return(col15)
    } else if (col_id == "col21") {
        return(col21)
    } else {
        stop(paste0("====> Error: wrong col_id"))
    }
}

#-- inputs --#
sample_prefix <- opt$prefix



trim_files <- vector()
merge_files <- vector()
stat_files <- vector()
canonical_barcode_files <- vector()
map_files  <- vector()
novel_barcode_files <- vector()



#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)



#-- reading files --#
total_reads <- vector()
merged_reads <- vector()
unmerged_reads <- vector()

inclusion_reads <- vector()
skipping_reads <- vector()
canonical_barcodes <- list()

map_reads <- vector()
unexplain_reads <- vector()
novel_barcodes <- list()


for(i in 1:length(reps))
{

}





#-- processing --#









#~~~~~~~~~~#
# function #
#~~~~~~~~~~#



#~~~~~~~~~~~~~~~#
# prepare files #
#~~~~~~~~~~~~~~~#
sample_prefix <- "HEK293T"
reps <- c(paste0(sample_prefix, "_Rep1"), paste0(sample_prefix, "_Rep2"), paste0(sample_prefix, "_Rep3"))

processing_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/1_processing"
trim_files <- vector()
merge_files <- vector()

filtering_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/2_filtering"
se_stat_files <- vector()
se_canonical_barcode_files <- vector()
pe_stat_files <- vector()
pe_canonical_barcode_files <- vector()

mapping_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/3_mapping"
se_map_files  <- vector()
se_novel_barcode_files <- vector()
pe_map_files  <- vector()
pe_novel_barcode_files <- vector()

for(i in 1:length(reps))
{
    trim_files <- append(trim_files, paste0(paste(processing_dir, reps[i], sep = "/"), "/", reps[i], ".trim.txt"))
    merge_files <- append(merge_files, paste0(paste(processing_dir, reps[i], sep = "/"), "/", reps[i], ".merge.txt"))

    se_stat_files <- append(se_stat_files, paste0(paste(filtering_dir, reps[i], sep = "/"), "/", reps[i], "_filter_se.filtered.idxstat.txt"))
    se_canonical_barcode_files <- append(se_canonical_barcode_files, paste0(paste(filtering_dir, reps[i], sep = "/"), "/", reps[i], "_filter_se.barcodes.txt"))
    pe_stat_files <- append(pe_stat_files, paste0(paste(filtering_dir, reps[i], sep = "/"), "/", reps[i], "_filter_pe.filtered.idxstat.txt"))
    pe_canonical_barcode_files <- append(pe_canonical_barcode_files, paste0(paste(filtering_dir, reps[i], sep = "/"), "/", reps[i], "_filter_pe.barcodes.txt"))

    se_map_files  <- append(se_map_files , paste0(paste(mapping_dir, reps[i], sep = "/"), "/", reps[i], "_map_se.summary.txt"))
    se_novel_barcode_files <- append(se_novel_barcode_files, paste0(paste(mapping_dir, reps[i], sep = "/"), "/", reps[i], "_map_se.barcodes.txt"))
    pe_map_files  <- append(pe_map_files , paste0(paste(mapping_dir, reps[i], sep = "/"), "/", reps[i], "_map_pe.summary.txt"))
    pe_novel_barcode_files <- append(pe_novel_barcode_files, paste0(paste(mapping_dir, reps[i], sep = "/"), "/", reps[i], "_map_pe.barcodes.txt"))
}

barcode_association_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/refseqs/"
barcode_association_file <- paste(barcode_association_dir, "varid_barcodes.unique.txt", sep = "/")

#~~~~~~~~~~~~#
# read files #
#~~~~~~~~~~~~#
total_reads <- vector()
merged_reads <- vector()
unmerged_reads <- vector()

se_inclusion_reads <- vector()
se_skipping_reads <- vector()
se_canonical_barcodes <- list()
pe_inclusion_reads <- vector()
pe_skipping_reads <- vector()
pe_canonical_barcodes <- list()

se_map_reads <- vector()
se_unexplain_reads <- vector()
se_novel_barcodes <- list()
pe_map_reads <- vector()
pe_unexplain_reads <- vector()
pe_novel_barcodes <- list()

for(i in 1:length(reps))
{
    total_reads <- append(total_reads, as.numeric(str_extract(readLines(trim_files[i], n = 2)[2], "\\d+")))
    merged_reads <- append(merged_reads, as.numeric(sub(".*Combined pairs:\\s*", "", grep("Combined pairs:", readLines(merge_files[i]), value = TRUE))))
    unmerged_reads <- append(unmerged_reads, as.numeric(sub(".*Uncombined pairs:\\s*", "", grep("Uncombined pairs:", readLines(merge_files[i]), value = TRUE))))

    se_inclusion_reads <- append(se_inclusion_reads, read.table(se_stat_files[i])[1, 3])
    se_skipping_reads <- append(se_skipping_reads, read.table(se_stat_files[i])[2, 3])
    se_canonical_barcodes[[reps[i]]] <- as.data.table(vroom(se_canonical_barcode_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    pe_inclusion_reads <- append(pe_inclusion_reads, read.table(pe_stat_files[i])[1, 3])
    pe_skipping_reads <- append(pe_skipping_reads, read.table(pe_stat_files[i])[2, 3])
    pe_canonical_barcodes[[reps[i]]] <- as.data.table(vroom(pe_canonical_barcode_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

    reads_tmp <- as.numeric(str_extract(readLines(se_map_files[i], n = 4)[4], "(?<=: )\\d+"))
    se_map_reads <- append(se_map_reads, reads_tmp)
    se_unexplain_reads <- append(se_unexplain_reads, as.numeric(str_extract(readLines(se_map_files[i], n = 2)[2], "(?<=: )\\d+")) - reads_tmp)
    reads_tmp <- as.numeric(str_extract(readLines(pe_map_files[i], n = 4)[4], "(?<=: )\\d+"))
    pe_map_reads <- append(pe_map_reads, reads_tmp)
    pe_unexplain_reads <- append(pe_unexplain_reads, as.numeric(str_extract(readLines(pe_map_files[i], n = 2)[2], "(?<=: )\\d+")) - reads_tmp)

    se_novel_barcodes[[reps[i]]] <- as.data.table(vroom(se_novel_barcode_files[i], delim = "\t", comment = "#", col_names = FALSE, show_col_types = FALSE))
    colnames(se_novel_barcodes[[reps[i]]]) <- c("var", "barcode")
    se_novel_barcodes[[reps[i]]] <- se_novel_barcodes[[reps[i]]][, .N, by = .(var, barcode)]
    setnames(se_novel_barcodes[[reps[i]]], "N", "count")
    pe_novel_barcodes[[reps[i]]] <- as.data.table(vroom(pe_novel_barcode_files[i], delim = "\t", comment = "#", col_names = FALSE, show_col_types = FALSE))
    colnames(pe_novel_barcodes[[reps[i]]]) <- c("var", "barcode")
    pe_novel_barcodes[[reps[i]]] <- pe_novel_barcodes[[reps[i]]][, .N, by = .(var, barcode)]
    setnames(pe_novel_barcodes[[reps[i]]], "N", "count")
}

barcode_association <- as.data.table(vroom(barcode_association_file, delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

#~~~~~~~~~#
# process #
#~~~~~~~~~#
summary_reads <- as.data.table(cbind(reps,
                                     total_reads, 
                                     merged_reads,
                                     unmerged_reads,
                                     se_inclusion_reads + pe_inclusion_reads,
                                     se_skipping_reads + pe_skipping_reads,
                                     se_map_reads + pe_map_reads,
                                     se_unexplain_reads + pe_unexplain_reads))
colnames(summary_reads) <- c("reps", 
                             "total_reads", 
                             "merged_reads", 
                             "unmerged_reads", 
                             "inclusion_reads", 
                             "skipping_reads",
                             "map_reads",
                             "unexplain_reads")
summary_reads_out <- paste0(sample_prefix, ".reads_stats.txt")
fwrite(summary_reads, summary_reads_out, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

summary_reads <- summary_reads %>%
                    mutate(across(c(total_reads, 
                                    merged_reads, 
                                    unmerged_reads, 
                                    inclusion_reads, 
                                    skipping_reads, 
                                    map_reads, 
                                    unexplain_reads), as.numeric))
summary_reads <- summary_reads %>%
                    mutate(merged_pct = round((merged_reads / total_reads) * 100, 1),
                           unmerged_pct = round((unmerged_reads / total_reads) * 100, 1),
                           library_coverage = as.integer(total_reads / length(unique(barcode_association$exon))),
                           inclusion_pct = round((inclusion_reads / total_reads) * 100, 1),
                           skipping_pct = round((skipping_reads / total_reads) * 100, 1),
                           map_pct = round((map_reads / total_reads) * 100, 1),
                           unexplain_pct = round((unexplain_reads / total_reads) * 100, 1))

tib_total_pct <- summary_reads %>% select(reps, merged_pct, unmerged_pct) %>%
                    pivot_longer(cols = -reps, names_to = "type", values_to = "pct")
tib_total_pct$type <- factor(tib_total_pct$type, levels = c("unmerged_pct", "merged_pct"))
tib_total_cov <- summary_reads %>% select(reps, library_coverage)
tib_total_cov$type <- "coverage"

select_colors <- select_colorblind("col15")[1:2]
fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)
y_scale <- max(tib_total_cov$library_coverage) * 2
png(paste0(sample_prefix, ".total_reads_pct.png"), width = 1200, height = 1600, units = "px", res = 300)
ggplot(tib_total_pct,  aes(x = reps, y = pct, fill = type)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_line(data = tib_total_cov, aes(x = reps, y = library_coverage / y_scale, group = 1), linetype = "dashed", color = "red", inherit.aes = FALSE) +
    geom_point(data = tib_total_cov, aes(x = reps, y = library_coverage / y_scale, color = type), shape = 18, size = 3, inherit.aes = FALSE) +
    scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~. * y_scale, name = "library coverage")) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = "red") +
    labs(x = NULL, y = "percent", title = sample_prefix) +
    theme(legend.position = "right", legend.title = element_blank()) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(aes(label = pct), position = position_fill(vjust = 0.5), size = 3)
dev.off()           
                        
tib_canonical_pct <- summary_reads %>% select(reps, inclusion_pct, skipping_pct) %>%
                        pivot_longer(cols = -reps, names_to = "type", values_to = "pct")

select_colors <- select_colorblind("col15")[1:2]
fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)
png(paste0(sample_prefix, ".canonical_reads_pct.png"), width = 1000, height = 1600, units = "px", res = 300)
ggplot(tib_canonical_pct,  aes(x = reps, y = pct, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors) +
    labs(x = NULL, y = "percent", title = sample_prefix) +
    theme(legend.position = "right", legend.title = element_blank()) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)
dev.off() 

tib_novel_pct <- summary_reads %>% select(reps, map_pct, unexplain_pct) %>%
                    pivot_longer(cols = -reps, names_to = "type", values_to = "pct")

select_colors <- select_colorblind("col15")[1:2]
fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)
png(paste0(sample_prefix, ".novel_reads_pct.png"), width = 1000, height = 1600, units = "px", res = 300)
ggplot(tib_novel_pct,  aes(x = reps, y = pct, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors) +
    labs(x = NULL, y = "percent", title = sample_prefix) +
    theme(legend.position = "right", legend.title = element_blank()) +
    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
    theme(axis.text = element_text(size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)
dev.off() 

barcode_can_ass_counts <- vector()
barcode_nol_ass_counts <- vector()
barcode_all_ass_counts <- vector()
variant_shared_counts <- vector()
for(i in 1:length(reps))
{
    venn_list <- list(unique(c(se_canonical_barcodes[[i]]$barcode, pe_canonical_barcodes[[i]]$barcode)), 
                      unique(c(se_novel_barcodes[[i]]$barcode, pe_novel_barcodes[[i]]$barcode)),
                      barcode_association$barcode)
    names(venn_list) <- c("canonical", "novel", "association")

    png(paste0(reps[i], ".barcodes_venn.png"), width = 1600, height = 1600, units = "px", res = 300)
    print(ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.2) + 
            scale_fill_gradient(low = "ivory", high = "tomato"))
    dev.off()

    venn_obj <- Venn(venn_list)
    venn_data <- process_data(venn_obj)
    barcode_can_ass_only <- venn_data$regionData %>% slice(5) %>% pull(item) %>% reduce(c)
    barcode_nol_ass_only <- venn_data$regionData %>% slice(6) %>% pull(item) %>% reduce(c)
    barcode_all_ass_only <- venn_data$regionData %>% slice(7) %>% pull(item) %>% reduce(c)

    barcode_can_ass_counts <- append(barcode_can_ass_counts, sum(se_canonical_barcodes[[i]][barcode %in% barcode_can_ass_only]$count) + 
                                                             sum(pe_canonical_barcodes[[i]][barcode %in% barcode_can_ass_only]$count))
    barcode_nol_ass_counts <- append(barcode_nol_ass_counts, sum(se_novel_barcodes[[i]][barcode %in% barcode_nol_ass_only]$count) +
                                                             sum(pe_novel_barcodes[[i]][barcode %in% barcode_nol_ass_only]$count))
    barcode_all_ass_counts <- append(barcode_all_ass_counts, sum(se_canonical_barcodes[[i]][barcode %in% barcode_all_ass_only]$count) + 
                                                             sum(pe_canonical_barcodes[[i]][barcode %in% barcode_all_ass_only]$count) +
                                                             sum(se_novel_barcodes[[i]][barcode %in% barcode_all_ass_only]$count) + 
                                                             sum(pe_novel_barcodes[[i]][barcode %in% barcode_all_ass_only]$count))

    barcode_shared <- unique(c(barcode_can_ass_only, barcode_nol_ass_only, barcode_all_ass_only))
    variant_shared_counts <- append(variant_shared_counts, length(unique(barcode_association[barcode %in% barcode_shared]$exon)))
}

barcode_counts <- as.data.table(cbind(reps,
                                      total_reads, 
                                      barcode_can_ass_counts,
                                      barcode_nol_ass_counts,
                                      barcode_all_ass_counts))
colnames(barcode_counts) <- c("reps", "total_reads", "canon_assoc", "novel_assoc", "canon_novel_assoc")
barcode_counts_out <- paste0(sample_prefix, ".barcodes_stats.txt")
fwrite(barcode_counts, barcode_counts_out, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

variant_counts <- as.data.table(cbind(reps,
                                      length(unique(barcode_association$exon)), 
                                      variant_shared_counts))
colnames(variant_counts) <- c("reps", "variants", "detected_variants")
variant_counts_out <- paste0(sample_prefix, ".variants_stats.txt")
fwrite(variant_counts, variant_counts_out, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

#~~~~~~~~#
# report #
#~~~~~~~~#
reads_image_list <- paste0("c(", "\"", paste0(sample_prefix, ".total_reads_pct.png"), "\",", 
                                 "\"", paste0(sample_prefix, ".canonical_reads_pct.png"), "\",", 
                                 "\"", paste0(sample_prefix, ".novel_reads_pct.png"), "\"", ")")

barcode_venn_list <- paste0("c(", "\"", paste0(reps[1], ".barcodes_venn.png"), "\",", 
                                  "\"", paste0(reps[2], ".barcodes_venn.png"), "\",", 
                                  "\"", paste0(reps[3], ".barcodes_venn.png"), "\"", ")")

junction_image_list <- paste0("c(", "\"", paste0(mapping_dir, "/", sample_prefix, ".junction_venn.png"), "\",", 
                                    "\"", paste0(mapping_dir, "/", sample_prefix, ".junction_corr.png"), "\",", 
                                    "\"", paste0(mapping_dir, "/", sample_prefix, ".junction_scatter.png"), "\"", ")")

report_path <- paste0(sample_prefix, ".splicing_report.Rmd")
sink(report_path)

cat("---", "\n", sep = "")
cat("title: \"Splicing Report\"", "\n", sep = "")
cat("date: \"`r format(Sys.time(), '%d %B %Y -- %A -- %X')`\"", "\n", sep = "")
cat("output:", "\n", sep = "")
cat("    html_document:", "\n", sep = "")
cat("        toc: true", "\n", sep = "")
cat("        toc_depth: 1", "\n", sep = "")
cat("        theme: united", "\n", sep = "")
cat("        highlight: tango", "\n", sep = "")
cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("```{r setup, include = FALSE}", "\n", sep = "")
cat("knitr::opts_chunk$set(echo = TRUE, fig.align = \"center\")", "\n", sep = "")
cat("library(reactable)", "\n", sep = "")
cat("```", "\n", sep = "")
cat("\n", sep = "")

cat("```{js, echo = FALSE}", "\n", sep = "")
cat("function formatNumber(num, precision = 1) {", "\n", sep = "")
cat("    const map = [", "\n", sep = "")
cat("        { suffix: 'T', threshold: 1e12 },", "\n", sep = "")
cat("        { suffix: 'B', threshold: 1e9 },", "\n", sep = "")
cat("        { suffix: 'M', threshold: 1e6 },", "\n", sep = "")
cat("        { suffix: 'K', threshold: 1e3 },", "\n", sep = "")
cat("        { suffix: '', threshold: 1 },", "\n", sep = "")
cat("    ];", "\n", sep = "")
cat("    const found = map.find((x) => Math.abs(num) >= x.threshold);", "\n", sep = "")
cat("    if (found) {", "\n", sep = "")
cat("        const formatted = (num / found.threshold).toFixed(precision) + found.suffix;", "\n", sep = "")
cat("        return formatted;", "\n", sep = "")
cat("    }", "\n", sep = "")
cat("    return num;", "\n", sep = "")
cat("}", "\n", sep = "")
cat("", "\n", sep = "")
cat("function rangeMore(column, state) {", "\n", sep = "")
cat("    let min = Infinity", "\n", sep = "")
cat("    let max = 0", "\n", sep = "")
cat("    state.data.forEach(function(row) {", "\n", sep = "")
cat("        const value = row[column.id]", "\n", sep = "")
cat("        if (value < min) {", "\n", sep = "")
cat("            min = Math.floor(value)", "\n", sep = "")
cat("        }", "\n", sep = "")
cat("        if (value > max) {", "\n", sep = "")
cat("            max = Math.ceil(value)", "\n", sep = "")
cat("        }", "\n", sep = "")
cat("    })", "\n", sep = "")
cat("", "\n", sep = "")
cat("    const filterValue = column.filterValue || min", "\n", sep = "")
cat("    const input = React.createElement('input', {", "\n", sep = "")
cat("        type: 'range',", "\n", sep = "")
cat("        value: filterValue,", "\n", sep = "")
cat("        min: min,", "\n", sep = "")
cat("        max: max,", "\n", sep = "")
cat("        onChange: function(event) {", "\n", sep = "")
cat("            column.setFilter(event.target.value || undefined)", "\n", sep = "")
cat("        },", "\n", sep = "")
cat("        style: { width: '100%', marginRight: '8px' },", "\n", sep = "")
cat("        'aria-label': 'Filter ' + column.name", "\n", sep = "")
cat("    })", "\n", sep = "")
cat("", "\n", sep = "")
cat("    return React.createElement(", "\n", sep = "")
cat("        'div',", "\n", sep = "")
cat("        { style: { display: 'flex', alignItems: 'center', height: '100%' } },", "\n", sep = "")
cat("        [input, formatNumber(filterValue)]", "\n", sep = "")
cat("    )", "\n", sep = "")
cat("}", "\n", sep = "")
cat("", "\n", sep = "")
cat("function filterMinValue(rows, columnId, filterValue) {", "\n", sep = "")
cat("    return rows.filter(function(row) {", "\n", sep = "")
cat("        return row.values[columnId] >= filterValue", "\n", sep = "")
cat("    })", "\n", sep = "")
cat("}", "\n", sep = "")
cat("", "\n", sep = "")
cat("function rangeLess(column, state) {", "\n", sep = "")
cat("    let min = Infinity", "\n", sep = "")
cat("    let max = 0", "\n", sep = "")
cat("    state.data.forEach(function(row) {", "\n", sep = "")
cat("        const value = row[column.id]", "\n", sep = "")
cat("        if (value < min) {", "\n", sep = "")
cat("            min = Math.floor(value)", "\n", sep = "")
cat("        }", "\n", sep = "")
cat("        if (value > max) {", "\n", sep = "")
cat("            max = Math.ceil(value)", "\n", sep = "")
cat("        }", "\n", sep = "")
cat("    })", "\n", sep = "")
cat("", "\n", sep = "")
cat("    const filterValue = column.filterValue || max", "\n", sep = "")
cat("    const input = React.createElement('input', {", "\n", sep = "")
cat("        type: 'range',", "\n", sep = "")
cat("        value: filterValue,", "\n", sep = "")
cat("        min: min,", "\n", sep = "")
cat("        max: max,", "\n", sep = "")
cat("        onChange: function(event) {", "\n", sep = "")
cat("            column.setFilter(event.target.value || undefined)", "\n", sep = "")
cat("        },", "\n", sep = "")
cat("        style: { width: '100%', marginRight: '8px' },", "\n", sep = "")
cat("        'aria-label': 'Filter ' + column.name", "\n", sep = "")
cat("    })", "\n", sep = "")
cat("", "\n", sep = "")
cat("    return React.createElement(", "\n", sep = "")
cat("        'div',", "\n", sep = "")
cat("        { style: { display: 'flex', alignItems: 'center', height: '100%' } },", "\n", sep = "")
cat("        [input, formatNumber(filterValue)]", "\n", sep = "")
cat("    )", "\n", sep = "")
cat("}", "\n", sep = "")
cat("", "\n", sep = "")
cat("function filterMaxValue(rows, columnId, filterValue) {", "\n", sep = "")
cat("    return rows.filter(function(row) {", "\n", sep = "")
cat("        return row.values[columnId] <= filterValue", "\n", sep = "")
cat("    })", "\n", sep = "")
cat("}", "\n", sep = "")
cat("```", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 1. Introduction", "\n", sep = "")
cat("This pipeline processes Illumina PE reads to detect the splicing events", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 2. Read Processing", "\n", sep = "")
cat("Displays plots and statistics for all samples.", "\n", sep = "")
cat("\n", sep = "")
cat("```{r, echo = FALSE}", "\n", sep = "")
cat("df <- as.data.frame(read.table(\"", summary_reads_out, "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150, align = \"left\"),", "\n", sep = "")
cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")))", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("```{r, echo = FALSE, fig.show = \"hold\", fig.align = \"default\", out.height = \"30%\", out.width = \"30%\"}", "\n", sep = "")
cat("knitr::include_graphics(", reads_image_list, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 3. Barcode Processing", "\n", sep = "")
cat("Displays plots and statistics for all samples.", "\n", sep = "")
cat("\n", sep = "")

cat("```{r, echo = FALSE, fig.show = \"hold\", fig.align = \"default\", out.height = \"32%\", out.width = \"32%\"}", "\n", sep = "")
cat("knitr::include_graphics(", barcode_venn_list, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 4. Junction Processing", "\n", sep = "")
cat("Displays plots and statistics for all samples.", "\n", sep = "")
cat("\n", sep = "")

cat("```{r, echo = FALSE, out.height = \"80%\", out.width = \"80%\"}", "\n", sep = "")
cat("knitr::include_graphics(", junction_image_list, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")

sink()
rmarkdown::render(report_path, clean = TRUE, quiet = TRUE)
invisible(file.remove(report_path))
