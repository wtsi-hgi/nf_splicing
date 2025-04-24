#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "UpSetR")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bed"),   type = "character", help = "input bed"),
                    make_option(c("-e", "--exon_pos"),    type = "character", help = "exon position file"),
                    make_option(c("-m", "--min_overlap"), type = "integer",   help = "min anchor for partial splicing", default = 3),
                    make_option(c("-o", "--output_dir"),  type = "character", help = "output directory",   default = getwd()))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

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
    for(i in 1:nrow(exon_pos))
    {
        if(start > exon_pos[i,]$exon_start & start < exon_pos[i,]$exon_end)
        {
            annotmp <- paste0("exon_splicing_3p_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        if(end > exon_pos[i,]$exon_start & end < exon_pos[i,]$exon_end)
        {
            annotmp <- paste0("exon_splicing_5p_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        if(start < exon_pos[i,]$exon_start & end > exon_pos[i,]$exon_end)
        {
            annotmp <- paste0("exon_skipping_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }
    }
    for(i in 1:nrow(intron_pos))
    {
        if(start > intron_pos[i,]$intron_start & start < intron_pos[i,]$intron_end)
        {
            annotmp <- paste0("intron_retension_5p_", intron_pos[i,]$intron_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        if(end > intron_pos[i,]$intron_start & end < intron_pos[i,]$intron_end)
        {
            annotmp <- paste0("intron_retension_3p_", intron_pos[i,]$intron_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }
    }

    return(annostr)
}









# initialization
exon_positions <- rbind(c("E10", 1, 210),
                        c("E11", 298, 444),
                        c("E12", 525, 690))
colnames(exon_positions) <- c("exon_id", "exon_start", "exon_end")
exon_positions <- as.data.table(exon_positions)
set(exon_positions, j = "exon_id", value = as.character(exon_positions$exon_id))
set(exon_positions, j = "exon_start", value = as.integer(exon_positions$exon_start))
set(exon_positions, j = "exon_end", value = as.integer(exon_positions$exon_end))

intron_positions <- rbind(c("I10", 211, 297),
                          c("I11", 445, 524))
colnames(intron_positions) <- c("intron_id", "intron_start", "intron_end")
intron_positions <- as.data.table(intron_positions)
set(intron_positions, j = "intron_id", value = as.character(intron_positions$intron_id))
set(intron_positions, j = "intron_start", value = as.integer(intron_positions$intron_start))
set(intron_positions, j = "intron_end", value = as.integer(intron_positions$intron_end))

categories <- c("intron_retention_I10",
                "intron_retention_I11",
                "intron_retension_5p_I10",
                "intron_retension_5p_I11",
                "intron_retension_3p_I10",
                "intron_retension_3p_I11",
                "exon_splicing_3p_E10",
                "exon_splicing_3p_E11",
                "exon_splicing_3p_E12",
                "exon_splicing_5p_E10",
                "exon_splicing_5p_E11",
                "exon_splicing_5p_E12",
                "exon_skipping_E11")

# read file
junction_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/3_mapping/HEK293T_Rep3"
junction_bed_file <- "HEK293T_Rep3.junctions.bed"
junction_cov_cutoff <- 2

junction_prefix <- strsplit(junction_bed_file, ".", fixed = T)[[1]][1]
setwd(junction_dir)

junction_bed <- as.data.table(read.table(junction_bed_file)) 

# convert bed to junction table
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

# classify junctions
annotate_junction <- function(start, end, exon_pos, intron_pos)
{
    if(start < min(exon_pos$exon_start) | end > max(exon_pos$exon_end))
    {
        return("out_of_range")
    }

    for(i in 1:nrow(intron_pos))
    {
        check_start <- (start >= (intron_pos[i,]$intron_start - 1) & start <= (intron_pos[i,]$intron_start + 1))
        check_end <- (end >= (intron_pos[i,]$intron_end - 1) & end <= (intron_pos[i,]$intron_end + 1))
        if(check_start & check_end)
        {
            if(i == 1)
            {
                return(paste0("intron_retention_", intron_pos[i + 1,]$intron_id))
            } else {
                return(paste0("intron_retention_", intron_pos[i - 1,]$intron_id))
            }
        }
    } 

    annostr <- ""
    annotmp <- ""
    for(i in 1:nrow(exon_pos))
    {
        if(start > exon_pos[i,]$exon_start & start < exon_pos[i,]$exon_end)
        {
            annotmp <- paste0("exon_splicing_3p_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        if(end > exon_pos[i,]$exon_start & end < exon_pos[i,]$exon_end)
        {
            annotmp <- paste0("exon_splicing_5p_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        if(start < exon_pos[i,]$exon_start & end > exon_pos[i,]$exon_end)
        {
            annotmp <- paste0("exon_skipping_", exon_pos[i,]$exon_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }
    }
    for(i in 1:nrow(intron_pos))
    {
        if(start > intron_pos[i,]$intron_start & start < intron_pos[i,]$intron_end)
        {
            annotmp <- paste0("intron_retension_5p_", intron_pos[i,]$intron_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }

        if(end > intron_pos[i,]$intron_start & end < intron_pos[i,]$intron_end)
        {
            annotmp <- paste0("intron_retension_3p_", intron_pos[i,]$intron_id)
            annostr <- ifelse(annostr == "", annotmp, paste0(annostr, ";", annotmp))
        }
    }

    return(annostr)
}

junction_annotation <- junction_data %>%
                       rowwise() %>%
                       mutate(Annotation = annotate_junction(Start, End, exon_positions, intron_positions)) %>%
                       ungroup()

write.table(junction_annotation, paste0(junction_prefix, ".junctions.txt"), sep = '\t', quote = FALSE, row.names = FALSE)

# junction_annotation <- read.table("I1_Rep1/I1_Rep1_map_se.junctions.txt", header = T, sep = "\t", stringsAsFactors = FALSE, quote = "")
# junction_annotation <- as_tibble(junction_annotation)
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

cov_filter <- function(row, cov_threshold)
{
    return(as.integer(row["Cov"]) > cov_threshold)
}

png(paste0(junction_prefix, ".junctions.png"), width = 1600, height = 1600, units = "px", res = 250)
upset(as.data.frame(junction_input), 
      nsets = length(categories), 
      order.by = "freq", 
      matrix.color = "yellowgreen", 
      main.bar.color = "royalblue", 
      sets.bar.color = "yellowgreen",
      queries = list(list(query = cov_filter, params = list(10), color = "tomato", active = TRUE)))
dev.off()

png(paste0(junction_prefix, ".variants.png"), width = 2000, height = 1600, units = "px", res = 250)
upset(as.data.frame(variant_input), 
      nsets = length(categories), 
      order.by = "freq", 
      matrix.color = "yellowgreen", 
      main.bar.color = "royalblue", 
      sets.bar.color = "yellowgreen",
      queries = list(list(query = cov_filter, params = list(10), color = "tomato", active = TRUE)))
dev.off()

# create unique junction table
junction_annotation_unique <- junction_annotation %>%
                              arrange(Start, End) %>%
                              mutate(group = cumsum(lag(Start, default = first(Start)) + 2 < Start | 
                                                    lag(End, default = first(End)) + 2 < End)) %>%
                              group_by(group) %>%
                              summarize(VarID = first(VarID),
                                        Cov = sum(as.numeric(Cov), na.rm = TRUE),
                                        Start = as.integer(median(Start, na.rm = TRUE)),
                                        End = as.integer(median(End, na.rm = TRUE)),
                                        Annotation = first(Annotation)) %>%
                              ungroup() %>%
                              select(-group) 
                           
write.table(junction_annotation_unique, paste0(junction_prefix, ".junctions.unique.txt"), sep = '\t', quote = FALSE, row.names = FALSE)
                           