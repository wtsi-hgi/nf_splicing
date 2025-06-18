#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "vroom", "gplots", "ggplot2", "ggVennDiagram", "scales", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-s", "--sample_id"),       type = "character", help = "list of sample IDs",                     default = NULL),
                    make_option(c("-n", "--novel_junctions"), type = "character", help = "list of novel classified junction file", default = NULL),
                    make_option(c("-e", "--exon_pos"),        type = "character", help = "exon position file",                     default = NULL),
                    make_option(c("-o", "--output_dir"),      type = "character", help = "output directory",                       default = getwd()),
                    make_option(c("-p", "--prefix"),          type = "character", help = "output prefix",                          default = NULL))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
  print_help(opt_parser)
  quit(status = 1)
}

# Check if required arguments are provided
if(is.null(opt$sample_id))       stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$novel_junctions)) stop("-n, list of novel junction files is required!", call. = FALSE)
if(is.null(opt$exon_pos))        stop("-e, exon position file is required!", call. = FALSE)

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
novel_junction_files <- unlist(strsplit(opt$novel_junctions, ","))

exon_positions <- read.table(opt$exon_pos, header = FALSE, sep = "\t")
colnames(exon_positions) <- c("exon_id", "exon_start", "exon_end")
exon_positions <- as.data.table(exon_positions)
set(exon_positions, j = "exon_id", value = as.character(exon_positions$exon_id))
set(exon_positions, j = "exon_start", value = as.integer(exon_positions$exon_start))
set(exon_positions, j = "exon_end", value = as.integer(exon_positions$exon_end))
exon_positions$length <- exon_positions$exon_end - exon_positions$exon_start + 1

intron_positions <- data.table(intron_id = paste0("I", 1:(nrow(exon_positions) - 1)),
                               intron_start = exon_positions$exon_end[1:(nrow(exon_positions) - 1)] + 1,
                               intron_end = exon_positions$exon_start[2:nrow(exon_positions)] - 1)
intron_positions[, length := intron_end - intron_start + 1]

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- opt$prefix

#-- reading files --#
novel_junctions <- list()
venn_junctions <- list()

for(i in 1:length(reps))
{
    novel_junctions[[reps[i]]] <- as.data.frame(vroom(novel_junction_files[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    novel_junctions[[reps[i]]] <- novel_junctions[[reps[i]]] %>% filter(Annotation != "out_of_range")
    
    venn_junctions[[reps[i]]] <- paste(novel_junctions[[reps[i]]]$VarID, novel_junctions[[reps[i]]]$Start, novel_junctions[[reps[i]]]$End, sep = "__")
}

venn_obj <- Venn(venn_junctions)
venn_data <- process_data(venn_obj)
junc_shared <- venn_data$regionData %>% slice(4:7) %>% pull(item) %>% reduce(c)

for(i in 1:length(reps))
{
    novel_junctions[[reps[i]]]$Index <- venn_junctions[[reps[i]]]
    novel_junctions[[reps[i]]] <- as_tibble(novel_junctions[[reps[i]]])
    novel_junctions[[reps[i]]] <- novel_junctions[[reps[i]]]
}

novel_junctions_shared <- lapply(novel_junctions, function(df) { df[df$Index %in% junc_shared, c("Index", "Cov")] })
for(i in 1:length(reps))
{
    if(i == 1)
    {
        join_tmp <- novel_junctions_shared[[i]]
    } else {
        join_tmp <- full_join(join_tmp, novel_junctions_shared[[i]], by = "Index")
    }
}

for(i in 1:length(reps))
{
    join_tmp <- join_tmp %>%
                    left_join(novel_junctions[[i]] %>% select(Index, Annotation), by = "Index") %>%
                    mutate(!!paste0("Annotation_R", i) := Annotation) %>%
                    select(-Annotation)
}

join_tmp <- join_tmp %>% 
                mutate(Merged_Col = coalesce(Annotation_R1, Annotation_R2, Annotation_R3)) %>%  
                select(-Annotation_R1, -Annotation_R2, -Annotation_R3)
colnames(join_tmp) <- c("Index", reps, "Annotation")
join_tmp <- join_tmp %>%
                separate(Index, into = c("VarID", "Start", "End"), sep = "__") %>%
                mutate(CovAvg = rowMeans(across(all_of(reps)), na.rm = TRUE))
join_tmp$Start <- as.numeric(join_tmp$Start)
join_tmp$End <- as.numeric(join_tmp$End)

write.table(join_tmp, paste0(sample_prefix, ".junction_category.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

cor_data <- join_tmp %>% select(all_of(reps))
cor_data[is.na(cor_data)] <- 0
cor_data_norm <- cor_data %>% mutate(across(everything(), ~ log2((. / sum(.)) * 1e6 + 1)))

#-- plotting --#
png(paste0(sample_prefix, ".junction_venn.png"), width = 1600, height = 1600, units = "px", res = 250)
ggVennDiagram(venn_junctions, label_alpha = 0, edge_size = 0.2) + 
    scale_fill_gradient(low = "ivory", high = "tomato")
dev.off()

png(paste0(sample_prefix, ".junction_scatter.png"), width = 1600, height = 1400, units = "px", res = 200)
ggplot(join_tmp, aes(x = Start, y = End, color = CovAvg)) +
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_vline(data = exon_positions, aes(xintercept = exon_start), color = "lightgrey", linetype = "dashed") +
    geom_vline(data = exon_positions, aes(xintercept = exon_end), color = "lightgrey", linetype = "dashed") +
    geom_hline(data = exon_positions, aes(yintercept = exon_start), color = "lightgrey", linetype = "dashed") +
    geom_hline(data = exon_positions, aes(yintercept = exon_end), color = "lightgrey", linetype = "dashed") +
    labs(title = sample_prefix, x = "Start", y = "End", size = "Average") +
    geom_point(shape = 19, alpha = 1, size = 0.5) +
    scale_color_continuous(trans = "log10", low = "blue", high = "brown1", labels = label_number(accuracy = 1)) +
    scale_x_continuous(limits = c(-10, max(exon_positions$exon_end)), breaks = seq(0, max(exon_positions$exon_end), by = 100)) +  
    scale_y_continuous(limits = c(-10, max(exon_positions$exon_end)), breaks = seq(0, max(exon_positions$exon_end), by = 100)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_rect(data = exon_positions, inherit.aes = FALSE, fill = "darkgreen", alpha = 0.6, aes(xmin = exon_start, xmax = exon_end, ymin = -10, ymax = 0)) + 
    geom_rect(data = exon_positions, inherit.aes = FALSE, fill = "darkgreen", alpha = 0.6, aes(ymin = exon_start, ymax = exon_end, xmin = -10, xmax = 0)) +
    geom_segment(data = intron_positions, inherit.aes = FALSE, fill = "grey", alpha = 0.6, aes(x = intron_start, xend = intron_end, y = -5, yend = -5)) +
    geom_segment(data = intron_positions, inherit.aes = FALSE, fill = "grey", alpha = 0.6, aes(y = intron_start, yend = intron_end, x = -5, xend = -5))
dev.off()

png(paste0(sample_prefix, ".junction_view.png"), width = 1600, height = 600, units = "px", res = 200)
ggplot() + 
    geom_rect(data = exon_positions, fill = t_col("darkgreen", 0.6), color = t_col("darkgreen", 0.6), aes(xmin = exon_start, xmax = exon_end, ymin = 0.8, ymax = 1)) +
    geom_rect(data = intron_positions, fill = "grey", color = "grey", aes(xmin = intron_start, xmax = intron_end, ymin = 0.87, ymax = 0.93)) +
    geom_curve(data = join_tmp, aes(x = Start, xend = End, y = 1, yend = 1, color = CovAvg, alpha = CovAvg), curvature = -0.5, linewidth = 0.2, lineend = "round") +
    scale_color_gradient(trans = "log10", low = "blue", high = "brown1", labels = label_number(accuracy = 1)) +
    scale_alpha(trans = "log10", range = c(0.01, 1)) + guides(alpha = "none") +
    coord_cartesian(xlim = c(0, max(exon_positions$exon_end)), ylim = c(0, 3)) +
    theme_void()
dev.off()

png(paste0(sample_prefix, ".junction_distribution.png"), width = 2000, height = 1200, units = "px", res = 200)
par(mfrow = c(2,1), mar = c(5, 4, 1, 1))
boxplot(join_tmp$CovAvg, log = 'x', horizontal = T, frame.plot = F, col = "forestgreen", outcol = "darkgrey", pch = 21, cex = 0.5, cex.axis = 0.8)
h <- hist(join_tmp$CovAvg, breaks = 40, plot = F)
h_log <- log10(h$counts + 1)
plot(h$mids, h_log, type = "h", lwd = 3, col = "forestgreen", xlab = "CovAvg", ylab = "Frequency", frame.plot = F, yaxt = "n", cex.axis = 0.8)
axis(side = 2, at = 0:max(ceiling(h_log)), labels = 10^(0:max(ceiling(h_log))), las = 1, cex.axis = 0.8)
dev.off()

png(paste0(sample_prefix, ".junction_corr.png"), width = 1600, height = 1600, units = "px", res = 100)
pairs(cor_data_norm,
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = function(x, y, ...) {panel.smooth(x, y, method = "lm", ...)},
      use = "complete.obs")
dev.off()
