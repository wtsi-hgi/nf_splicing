#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "vroom", "ggVennDiagram", "htmltools", "reactable", "optparse", "sparkline", "UpSetR")
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
                    make_option(c("-j", "--junction_plots"),      type = "character", help = "list of junction plots",                    default = NULL),
                    make_option(c("-g", "--junction_category"),   type = "character", help = "junction category file",                    default = NULL),
                    make_option(c("-d", "--splicing_matrices"),   type = "character", help = "list of splicing matrices",                 default = NULL),
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
if(is.null(opt$junction_plots))      stop("-j, list of junction plots is required!", call. = FALSE)
if(is.null(opt$junction_category))   stop("-g, junction category file is required!", call. = FALSE)
if(is.null(opt$splicing_matrices))   stop("-d, list of splicing matrices is required!", call. = FALSE)

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


panel.hist <- function(x, ...)
{
    x <- na.omit(x)

    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    his <- hist(x, plot = FALSE)
    breaks <- his$breaks
    nB <- length(breaks)
    y <- his$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = t_col("yellowgreen", 0.8), ...)
    lines(density(x), col = 2, lwd = 2)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))

    Cor <- abs(cor(x, y, use = "complete.obs"))
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])

    if(missing(cex.cor))
    {
        cex.cor <- 0.4 / strwidth(txt)
    }

    # Generate heatmap color scale based on correlation value
    col.range = c("blue", "white", "red")
    heatmap_col <- colorRampPalette(col.range)(100)
    col_index <- round(Cor * 100)
    bg.col <- heatmap_col[col_index]

    # Draw the rectangle with heatmap color
    rect(0, 0, 1, 1, col = bg.col, border = NA)
    
    # Draw correlation text on top
    text(0.5, 0.5, txt, cex = 0.6 + cex.cor * Cor)
}

panel.smooth <- function(x, y, 
                         col.line = "red", 
                         plot.type = c("point", "smooth"),
                         method = c("loess", "lm"), span = 0.5, degree = 2, level = 0.95, 
                         pch = 21, col.pch = t_col("royalblue", 0.8), col.bg = t_col("royalblue", 0.4), cex.pch = 0.8, ...)
{
    plot.type <- match.arg(plot.type, c("point", "smooth"))

    data <- na.omit(data.frame(x, y))
    x <- data$x
    y <- data$y

    if(plot.type == "point")
    {
        points(x, y, pch = pch, col = col.pch, bg = col.bg, cex = cex.pch, ...)
    } else {
        smoothScatter(x, y, add = TRUE, nrpoints = 0)
    }

    method <- match.arg(method)
    data <- subset(data, x != 0 & y != 0)
    x <- data$x
    y <- data$y
    if(method == "loess")
    {
        loess_fit <- loess(y ~ x, span = span, degree = degree)
        pred <- predict(loess_fit, se = TRUE)
        lines(x, pred$fit, col = col.line, ...)
    } else if (method == "lm") {
        lm_fit <- lm(y ~ x)
        pred <- predict(lm_fit, interval = "confidence", level = level)
        lines(x, pred[, "fit"], col = col.line, ...)
    }
}

#-- inputs --#
reps <- unlist(strsplit(opt$sample_id, ","))
trim_files <- unlist(strsplit(opt$trim_stats, ","))
merge_files <- unlist(strsplit(opt$merge_stats, ","))
filter_files <- unlist(strsplit(opt$filter_idxstats, ","))
map_files  <- unlist(strsplit(opt$map_stats, ","))

canonical_barcode_files <- unlist(strsplit(opt$canonical_barcodes, ","))
novel_barcode_files <- unlist(strsplit(opt$novel_barcodes, ","))

splicing_matrices_files <- unlist(strsplit(opt$splicing_matrices, ","))

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- opt$prefix

#-- reading files --#
barcode_association <- as.data.table(vroom(opt$barcode_association, delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

total_reads <- vector()
merged_reads <- vector()
unmerged_reads <- vector()

inclusion_reads <- vector()
skipping_reads <- vector()
canonical_barcodes <- list()

map_reads <- vector()
unexplain_reads <- vector()
novel_barcodes <- list()

splicing_matrices <- list()

for(i in 1:3)
{
    tmp_value <- as.numeric(str_extract(grep("reads passed filter:", readLines(trim_files[i]), value = TRUE), "\\d+"))
    total_reads <- append(total_reads, tmp_value / 2)

    tmp_value <- as.numeric(str_extract(grep("Combined pairs:", readLines(merge_files[i]), value = TRUE), "\\d+"))
    merged_reads <- append(merged_reads, tmp_value)

    tmp_value <- as.numeric(str_extract(grep("Uncombined pairs:", readLines(merge_files[i]), value = TRUE), "\\d+"))
    unmerged_reads <- append(unmerged_reads, tmp_value)

    inclusion_reads <- append(inclusion_reads, as.numeric(read.table(filter_files[i])[1, 2]))
    skipping_reads <- append(skipping_reads, as.numeric(read.table(filter_files[i])[2, 2]))
    canonical_barcodes[[reps[i]]] <- as.data.table(vroom(canonical_barcode_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

    map_reads <- append(map_reads, as.numeric(read.table(map_files[i])[1, 2]))
    unexplain_reads <- append(unexplain_reads, as.numeric(read.table(map_files[i])[2, 2]) - as.numeric(read.table(map_files[i])[1, 2]))
    novel_barcodes[[reps[i]]] <- as.data.table(vroom(novel_barcode_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))

    splicing_matrices[[reps[i]]] <- as.data.table(vroom(splicing_matrices_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
}

#-- processing --#
summary_reads <- as.data.table(cbind(reps,
                                     total_reads, 
                                     merged_reads,
                                     unmerged_reads,
                                     inclusion_reads,
                                     skipping_reads,
                                     map_reads,
                                     unexplain_reads))
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
                           library_coverage = as.integer(total_reads / length(unique(barcode_association$varid))),
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
png(paste0(sample_prefix, ".total_reads_pct.png"), width = 1100, height = 1600, units = "px", res = 300)
ggplot(tib_total_pct,  aes(x = reps, y = pct, fill = type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
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
    venn_list <- list(unique(canonical_barcodes[[i]]$barcode), 
                      unique(novel_barcodes[[i]]$barcode),
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

    barcode_can_ass_counts <- append(barcode_can_ass_counts, sum(canonical_barcodes[[i]][barcode %in% barcode_can_ass_only]$count))
    barcode_nol_ass_counts <- append(barcode_nol_ass_counts, sum(novel_barcodes[[i]][barcode %in% barcode_nol_ass_only]$count))
    barcode_all_ass_counts <- append(barcode_all_ass_counts, sum(canonical_barcodes[[i]][barcode %in% barcode_all_ass_only]$count) + 
                                                             sum(novel_barcodes[[i]][barcode %in% barcode_all_ass_only]$count))

    barcode_shared <- unique(c(barcode_can_ass_only, barcode_nol_ass_only, barcode_all_ass_only))
    variant_shared_counts <- append(variant_shared_counts, length(unique(barcode_association[barcode %in% barcode_shared]$varid)))
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
                                      length(unique(barcode_association$varid)), 
                                      variant_shared_counts))
colnames(variant_counts) <- c("reps", "variants", "detected_variants")
variant_counts_out <- paste0(sample_prefix, ".variants_stats.txt")
fwrite(variant_counts, variant_counts_out, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

junction_category <- as_tibble(vroom(opt$junction_category, delim = "\t", col_names = TRUE, show_col_types = FALSE))
junction_category_reshape <- junction_category %>%
                                mutate(CovAvg_log = log2(CovAvg + 1)) %>%
                                separate_rows(Annotation, sep = ";") %>%
                                select(c(VarID, Annotation, CovAvg_log))
upset_input <- junction_category_reshape %>%
                distinct(VarID, Annotation) %>%
                mutate(value = 1) %>%
                pivot_wider(names_from = Annotation, values_from = value, values_fill = 0)
upset_join <- junction_category_reshape %>% left_join(upset_input, by = "VarID")
png(paste0(sample_prefix, ".junction_category.png"), width = 2400, height = 1600, units = "px", res = 200)
upset(as.data.frame(upset_join), 
      sets = unique(upset_join$Annotation), 
      order.by = "freq", 
      matrix.color = "yellowgreen",
      main.bar.color = "royalblue", 
      sets.bar.color = "yellowgreen",
      boxplot.summary = c("CovAvg_log"))
dev.off()

for(i in 1:length(reps))
{
    splicing_matrices[[reps[i]]] <- splicing_matrices[[reps[i]]] %>%
                                        rowwise() %>%
                                        mutate(PSI = canonical_inclusion_E2 / sum(c_across(-varid))) %>%
                                        ungroup()
}

psi_list <- lapply(names(splicing_matrices), function(rep_id) { splicing_matrices[[rep_id]] %>%
                                                                select(varid, PSI) %>%
                                                                rename(!!rep_id := PSI) })
psi_wide <- reduce(psi_list, full_join, by = "varid")

png(paste0(sample_prefix, ".psi_correlation.png"), width = 1200, height = 1200, units = "px", res = 100)
pairs(psi_wide[, -1],
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = function(x, y, ...) {panel.smooth(x, y, method = "lm", ...)},
      use = "complete.obs")
dev.off()

#-- reporting --#
reads_image_list <- paste0("c(", "\"", paste0(sample_prefix, ".total_reads_pct.png"), "\",", 
                                 "\"", paste0(sample_prefix, ".canonical_reads_pct.png"), "\",", 
                                 "\"", paste0(sample_prefix, ".novel_reads_pct.png"), "\"", ")")

barcode_venn_list <- paste0("c(", "\"", paste0(reps[1], ".barcodes_venn.png"), "\",", 
                                  "\"", paste0(reps[2], ".barcodes_venn.png"), "\",", 
                                  "\"", paste0(reps[3], ".barcodes_venn.png"), "\"", ")")

junction_plots <- strsplit(opt$junction_plots, ",")[[1]]
junction_image_list1 <- "c("
for (i in 1:2) {
    junction_image_list1 <- paste0(junction_image_list1, "\"", junction_plots[i], "\"")
    if (i < 2) {
        junction_image_list1 <- paste0(junction_image_list1, ",")
    } else {
        junction_image_list1 <- paste0(junction_image_list1, ")")
    }
}

junction_image_list2 <- "c("
for (i in 3:length(junction_plots)) {
    junction_image_list2 <- paste0(junction_image_list2, "\"", junction_plots[i], "\"")
    if (i < length(junction_plots)) {
        junction_image_list2 <- paste0(junction_image_list2, ",")
    } else {
        junction_image_list2 <- paste0(junction_image_list2, ")")
    }
}

report_path <- paste0(sample_prefix, ".splicing_report.Rmd")
sink(report_path)

cat("---", "\n", sep = "")
cat("title: \"Splicing Report\"", "\n", sep = "")
cat("date: \"`r format(Sys.time(), '%d %B %Y -- %A -- %X')`\"", "\n", sep = "")
cat("output:", "\n", sep = "")
cat("    html_document:", "\n", sep = "")
cat("        toc: true", "\n", sep = "")
cat("        toc_depth: 4", "\n", sep = "")
cat("        theme: united", "\n", sep = "")
cat("        highlight: tango", "\n", sep = "")
cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("```{r setup, include = FALSE}", "\n", sep = "")
cat("knitr::opts_chunk$set(echo = TRUE, fig.align = \"center\")", "\n", sep = "")
cat("library(reactable)", "\n", sep = "")
cat("library(sparkline)", "\n", sep = "")
cat("library(UpSetR)", "\n", sep = "")
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
cat("This pipeline is designed to quantify splicing events within minigene-based mutagenesis libraries. 
     Minigene systems are widely used to study the regulatory mechanisms of RNA splicing, 
     particularly in the context of variant interpretation and functional genomics. 
     By introducing specific mutations into synthetic constructs (minigenes), 
     researchers can assess how sequence changes affect splicing outcomes in a controlled cellular environment.", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 2. Read Processing", "\n", sep = "")
cat("This section summarises the distribution of reads according to the alignments.", "\n", sep = "")
cat("\n", sep = "")
cat("* **total_reads:** the number of raw reads aftering adaptor trimming and quality filtering.", "\n", sep = "")
cat("* **merged_reads:** the number of reads supporting the splicing events.", "\n", sep = "")
cat("* **unmerged_reads:** the number of reads supporting the whole minigene transcript without any splicing event.", "\n", sep = "")
cat("* **inclusion_reads:** the number of reads supporting canonical exon inclusion.", "\n", sep = "")
cat("* **skipping_reads:** the number of reads supporting canonical exon skipping.", "\n", sep = "")
cat("* **map_reads:** the number of reads supporting novel splicing events.", "\n", sep = "")
cat("* **unexplain_reads:** the number of reads which cannot be categorized into any event.", "\n", sep = "")
cat("\n", sep = "")
cat("```{r, echo = FALSE}", "\n", sep = "")
cat("df <- as.data.frame(read.table(\"", summary_reads_out, "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
cat("min_row <- ifelse(nrow(df) > 10, 10, nrow(df))", "\n", sep = "")
cat("reactable(df, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
cat("          filterable = TRUE, minRows = min_row, defaultColDef = colDef(minWidth = 150, align = \"left\"))", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("```{r, echo = FALSE, fig.show = \"hold\", fig.align = \"default\", out.height = \"30%\", out.width = \"30%\"}", "\n", sep = "")
cat("knitr::include_graphics(", reads_image_list, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 3. Barcode Processing", "\n", sep = "")
cat("Pipeline include a method of barcode detection which may differ from the method of input barcode association.
     This section shows the discrepancy between the pipeline barcodes and the input barcodes", "\n", sep = "")
cat("\n", sep = "")

cat("```{r, echo = FALSE, fig.show = \"hold\", fig.align = \"default\", out.height = \"32%\", out.width = \"32%\"}", "\n", sep = "")
cat("knitr::include_graphics(", barcode_venn_list, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("## 4. Junction Processing", "\n", sep = "")
cat("This section summarises all the novel splicing events", "\n", sep = "")
cat("\n", sep = "")

cat("### 4.1. Correlations between replicates", "\n", sep = "")
cat("\n", sep = "")
cat("```{r, echo = FALSE, fig.show = \"hold\", fig.align = \"default\", out.height = \"48%\", out.width = \"48%\"}", "\n", sep = "")
cat("knitr::include_graphics(", junction_image_list1, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("### 4.2. Distribution of splicing junctions (appearing in 2 replicates at least)", "\n", sep = "")
cat("\n", sep = "")
cat("```{r, echo = FALSE, out.height = \"80%\", out.width = \"80%\"}", "\n", sep = "")
cat("knitr::include_graphics(", junction_image_list2, ", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("---", "\n", sep = "")
cat("\n", sep = "")

cat("### 4.3. Splicing events  (appearing in 2 replicates at least) by variants", "\n", sep = "")
cat("\n", sep = "")
cat("```{r, echo = FALSE}", "\n", sep = "")
cat("junction_category <- as_tibble(vroom(\"", opt$junction_category, "\", delim = \"\\t\", col_names = TRUE, show_col_types = FALSE))", "\n", sep = "")
cat("junction_category_summary <- junction_category %>%", "\n", sep = "")
cat("                               separate_rows(Annotation, sep = \";\") %>%", "\n", sep = "")
cat("                               group_by(Annotation) %>%", "\n", sep = "")
cat("                               summarise(CovAvg_log = list(log10(CovAvg + 1)),", "\n", sep = "")
cat("                                         No_of_variant = n_distinct(VarID),", "\n", sep = "")
cat("                                         .groups = \"drop\")", "\n", sep = "")
cat("boxplot_min <- min(unlist(junction_category_summary$CovAvg_log), na.rm = TRUE)", "\n", sep = "")
cat("boxplot_max <- max(unlist(junction_category_summary$CovAvg_log), na.rm = TRUE)", "\n", sep = "")
cat("reactable(junction_category_summary, highlight = TRUE, bordered = TRUE, striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
cat("          filterable = FALSE, sortable = FALSE, minRows = 12, defaultColDef = colDef(minWidth = 200, align = \"left\"), defaultPageSize = 12,", "\n", sep = "")
cat("          columns = list(", "\n", sep = "")
cat("              Annotation = colDef(name = \"Splicing Type\"),", "\n", sep = "")
cat("              CovAvg_log = colDef(name = \"log10(CovAvg+1)\",", "\n", sep = "")
cat("                                  cell = function(value, index) {", "\n", sep = "")
cat("                                             if (length(value) > 5) {", "\n", sep = "")
cat("                                                 sparkline(value, type = \"box\", width = 180,", "\n", sep = "")
cat("                                                           chartRangeMin = boxplot_min, chartRangeMax = boxplot_max) }}),", "\n", sep = "")
cat("              No_of_variant = colDef(name = \"No. of Variants\")))", "\n", sep = "")
cat("", "\n", sep = "")
cat("knitr::include_graphics(\"", paste0(sample_prefix, ".junction_category.png"), "\", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")
cat("<br>", "\n", sep = "")
cat("\n", sep = "")

cat("## 5. PSI Correlation", "\n", sep = "")
cat("This section summarises the correlation of PSI values between replicates.", "\n", sep = "")

cat("\n", sep = "")
cat("```{r, echo = FALSE, fig.align = \"center\", out.width = \"50%\"}", "\n", sep = "")
cat("knitr::include_graphics(\"", paste0(sample_prefix, ".psi_correlation.png"), "\", rel_path = FALSE)", "\n", sep = "")
cat("```", "\n", sep = "")

sink()
rmarkdown::render(report_path, clean = TRUE, quiet = TRUE)
invisible(file.remove(report_path))
