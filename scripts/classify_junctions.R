#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "UpSetR", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bed"),   type = "character", help = "input bed",                        default = NULL),
                    make_option(c("-e", "--exon_pos"),    type = "character", help = "exon position file",               default = NULL),
                    make_option(c("-m", "--min_overlap"), type = "integer",   help = "min anchor for partial splicing",  default = 2),
                    make_option(c("-c", "--min_cov"),     type = "integer",   help = "minimum coverage cutoff",          default = 2),
                    make_option(c("-r", "--reduce"),      type = "integer",   help = "the number of bases for reducing", default = 2),
                    make_option(c("-o", "--output_dir"),  type = "character", help = "output directory",                 default = getwd()),
                    make_option(c("-p", "--prefix"),      type = "character", help = "output prefix",                    default = NULL))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
  print_help(opt_parser)
  quit(status = 1)
}

# Check if required arguments are provided
if(is.null(opt$input_bed)) stop("-b, input bed file is required!", call. = FALSE)
if(is.null(opt$exon_pos))  stop("-e, exon position file is required!", call. = FALSE)

#-- function --#
annotate_junction <- function(start, end, exon_pos, intron_pos)
{
    if(start < min(exon_pos$exon_start) | end > max(exon_pos$exon_end)) return("out_of_range")

    for(i in 1:nrow(intron_pos))
    {
        # intron retention: allow 1bp mismatch for start and end
        check_start <- (start >= (intron_pos[i,]$intron_start - 1) & start <= (intron_pos[i,]$intron_start + 1))
        check_end <- (end >= (intron_pos[i,]$intron_end - 1) & end <= (intron_pos[i,]$intron_end + 1))
        if(check_start & check_end)
        {
            ifelse(i == 1, 
                   return(paste0("intron_retention_", intron_pos[i + 1,]$intron_id)), 
                   return(paste0("intron_retention_", intron_pos[i - 1,]$intron_id)))
        }
    } 

    annostr <- ""
    annotmp <- ""
    # for exon splicing
    for(i in 1:nrow(exon_pos))
    {
        cond1 <- start > (exon_pos[i,]$exon_start + opt$min_overlap)
        cond2 <- start < (exon_pos[i,]$exon_end - opt$min_overlap)
        if(cond1 & cond2)
        {
            annotmp <- paste0("exon_splicing_3p_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        cond1 <- end > (exon_pos[i,]$exon_start + opt$min_overlap)
        cond2 <- end < (exon_pos[i,]$exon_end - opt$min_overlap)
        if(cond1 & cond2)
        {
            annotmp <- paste0("exon_splicing_5p_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        cond1 <- start < (exon_pos[i,]$exon_start + opt$min_overlap)
        cond2 <- end > (exon_pos[i,]$exon_end - opt$min_overlap)
        if(cond1 & cond2)
        {
            annotmp <- paste0("exon_skipping_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }
    }
    # for intron splicing
    for(i in 1:nrow(intron_pos))
    {
        cond1 <- start > intron_pos[i,]$intron_start + opt$min_overlap
        cond2 <- start < intron_pos[i,]$intron_end - opt$min_overlap
        if(cond1 & cond2)
        {
            annotmp <- paste0("intron_retension_5p_", intron_pos[i,]$intron_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        cond1 <- end > intron_pos[i,]$intron_start + opt$min_overlap
        cond2 <- end < intron_pos[i,]$intron_end - opt$min_overlap
        if(cond1 & cond2)
        {
            annotmp <- paste0("intron_retension_3p_", intron_pos[i,]$intron_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }
    }

    return(annostr)
}

cov_filter <- function(row, cov_threshold)
{
    return(as.integer(row["Cov"]) > cov_threshold)
}

#-- inputs --#
junction_bed_file <- opt$input_bed
junction_cov_cutoff <- opt$min_cov

junction_prefix <- ifelse(is.null(opt$prefix), 
                          tools::file_path_sans_ext(basename(junction_bed_file)), 
                          opt$prefix)

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

categories <- c(paste0("intron_retention_", intron_positions$intron_id),
                paste0("intron_retension_5p_", intron_positions$intron_id),
                paste0("intron_retension_3p_", intron_positions$intron_id),
                paste0("exon_splicing_3p_", exon_positions$exon_id),
                paste0("exon_splicing_5p_", exon_positions$exon_id),
                paste0("exon_skipping_", exon_positions$exon_id[2]))

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

output_junction <- paste0(junction_prefix, ".classified_junctions.txt")
if(file.exists(output_junction)) invisible(file.remove(output_junction))

png_junction <- paste0(junction_prefix, ".classified_junctions.png")
if(file.exists(png_junction)) invisible(file.remove(png_junction))

png_variant <- paste0(junction_prefix, ".classified_variants.png")
if(file.exists(png_variant)) invisible(file.remove(png_variant))

output_junction_reduce <- paste0(junction_prefix, ".classified_junctions.reduce.txt")
if(file.exists(output_junction_reduce)) invisible(file.remove(output_junction_reduce))

#-- processing --#
junction_bed <- fread(junction_bed_file, header = FALSE, sep = '\t')

junction_data <- junction_bed %>%
                 separate(V11, into = c("block_length1", "block_length2"), sep = ",") %>%
                 separate(V12, into = c("block_start1", "block_start2"), sep = ",") %>%
                 mutate(across(starts_with("block_"), as.numeric)) %>%
                 mutate(JuncStart = as.numeric(V2) + block_start1 + block_length1 + 1,
                        JuncEnd = as.numeric(V2) + block_start2 - 1) %>%
                 select(Column1 = V1, Column5 = V5, JuncStart, JuncEnd) %>%
                 rename(VarID = Column1, Cov = Column5, Start = JuncStart, End = JuncEnd) %>%
                 group_by(VarID, Start, End) %>%
                 summarise(Cov = sum(Cov), .groups = "drop") %>%
                 filter(as.numeric(Cov) >= junction_cov_cutoff) %>%
                 as_tibble()

junction_annotation <- junction_data %>%
                       rowwise() %>%
                       mutate(Annotation = annotate_junction(Start, End, exon_positions, intron_positions)) %>%
                       ungroup()

write.table(junction_annotation, output_junction, sep = '\t', quote = FALSE, row.names = FALSE)

junction_input <- junction_annotation %>%
                  separate_rows(Annotation, sep = ";", convert = TRUE) %>%
                  mutate(Annotation = str_trim(Annotation)) %>%
                  mutate(value = 1) %>%
                  distinct() %>%  # remove duplicated rows
                  pivot_wider(names_from = Annotation, values_from = value, values_fill = 0) %>%
                  { if ("out_of_range" %in% colnames(.)) select(., -out_of_range) else . }
missing_categories <- setdiff(categories, colnames(junction_input))
junction_input[missing_categories] <- 0
junction_input <- junction_input[, c(colnames(junction_input)[1:4], categories)]

variant_input <- junction_input %>%
                 group_by(VarID) %>%
                 summarise(Cov = sum(Cov, na.rm = TRUE),
                           across(all_of(categories), ~ as.integer(any(.x == 1, na.rm = TRUE)))) %>%
                 ungroup()

png(png_junction, width = 1600, height = 1600, units = "px", res = 250)
plot(1:10, 1:10, type = "n", xlab = "", ylab = "")
# upset(as.data.frame(junction_input), 
#       nsets = length(categories), 
#       order.by = "freq", 
#       matrix.color = "yellowgreen", 
#       main.bar.color = "royalblue", 
#       sets.bar.color = "yellowgreen",
#       queries = list(list(query = cov_filter, params = list(10), color = "tomato", active = TRUE)))
dev.off()

png(png_variant, width = 2000, height = 1600, units = "px", res = 250)
plot(1:10, 1:10, type = "n", xlab = "", ylab = "")
# upset(as.data.frame(variant_input), 
#       nsets = length(categories), 
#       order.by = "freq", 
#       matrix.color = "yellowgreen", 
#       main.bar.color = "royalblue", 
#       sets.bar.color = "yellowgreen",
#       queries = list(list(query = cov_filter, params = list(10), color = "tomato", active = TRUE)))
dev.off()

# merge junctions with similar start and end positions to reduce the number of junctions
# eg: Junc1: 100-200, Junc2: 101-201, Junc3: 102-202
# merge them into one junction: 100-202
junction_annotation_reduce<- junction_annotation %>%
                              arrange(Start, End) %>%
                              mutate(group = cumsum(lag(Start, default = first(Start)) + opt$reduce < Start | 
                                                    lag(End, default = first(End)) + opt$reduce < End)) %>%
                              group_by(group) %>%
                              summarize(VarID = first(VarID),
                                        Cov = sum(as.numeric(Cov), na.rm = TRUE),
                                        Start = as.integer(median(Start, na.rm = TRUE)),
                                        End = as.integer(median(End, na.rm = TRUE)),
                                        Annotation = first(Annotation)) %>%
                              ungroup() %>%
                              select(-group) 
                           
write.table(junction_annotation_reduce, output_junction_reduce, sep = '\t', quote = FALSE, row.names = FALSE)
                           