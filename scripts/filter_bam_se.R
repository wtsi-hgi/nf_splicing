#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Rsamtools", "Biostrings", "parallel", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bam"),  type = "character", help = "input bam",          default = NULL),
                    make_option(c("-e", "--exon_pos"),   type = "character", help = "exon position file", default = NULL),
                    make_option(c("-s", "--softclip"),   type = "integer",   help = "softclip tolerance", default = 5),
                    make_option(c("-o", "--output_dir"), type = "character", help = "output directory",   default = getwd()),
                    make_option(c("-t", "--threads"),    type = "integer",   help = "number of threads",  default = 64),
                    make_option(c("-c", "--chunk_size"), type = "integer",   help = "chunk size",         default = 10000))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if(is.null(opt$input_bam)) stop("-b, input bam file is required!", call. = FALSE)
if(is.null(opt$exon_pos))  stop("-e, exon position file is required!", call. = FALSE)

#-- function --#
calc_end_pos <- function(start_pos, cigar)
{
    tb_cigar <- tibble(cigar = cigar) %>%
                  mutate(values = str_extract_all(cigar, "\\d+"), operations = str_extract_all(cigar, "[A-Z]")) %>%  
                  unnest(c(values, operations)) %>%  # Expand the list columns
                  mutate(value = as.numeric(values)) %>%
                  select(-values)

    valid_ops <- c("M", "D", "N")

    total_ref_length <- tb_cigar %>%
                          filter(operations %in% valid_ops) %>%
                          summarise(ref_length = sum(value, na.rm = TRUE)) %>%
                          pull(ref_length)
    
    end_pos <- start_pos + total_ref_length - 1
    return(end_pos)
}

write_to_sam_file <- function(bam_data, index, output_sam)
{
    read_out <- c(bam_data$qname[index],
                  bam_data$flag[index],
                  as.character(bam_data$rname[index]),
                  bam_data$pos[index],
                  bam_data$mapq[index],
                  bam_data$cigar[index],
                  ifelse(is.na(bam_data$mrnm[index]), "*", as.character(bam_data$mrnm[index])),
                  ifelse(is.na(bam_data$mpos[index]), 0, bam_data$mpos[index]),
                  bam_data$isize[index],
                  data.frame(bam_data$seq[index])[1,1],
                  data.frame(bam_data$qual[index])[1,1],
                  paste0("NM:i:", bam_data$tag$NM[index]),
                  paste0("AS:i:", bam_data$tag$AS[index]),
                  paste0("XS:i:", bam_data$tag$XS[index]))

    fwrite(as.data.frame(t(read_out)), output_sam, append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

process_read <- function(i, bam_data)
{
    read_name <- bam_data[[1]]$qname[i]
    read_start_pos <- bam_data[[1]]$pos[i]
    read_cigar <- bam_data[[1]]$cigar[i]
    read_end_pos <- calc_end_pos(read_start_pos, read_cigar)
    read_ref <- bam_data[[1]]$rname[i]
    read_seq <- data.frame(bam_data[[1]]$seq[i])[1,1]
    read_qual <- data.frame(bam_data[[1]]$qual[i])[1,1]

    # Note:
    # 1. read alignment should be end-to-end, allow 5 bp difference for softclip
    # 2. library may be longer than reads, so reads may not cover the end part of last exon, allow 5 bp difference
    read_start_pass <- ((exon_positions$exon_end[1] - read_start_pos) > pos_tolerance)
    read_end_pass <- ifelse(read_ref == "E1_E2_E3", 
                            ((read_end_pos - sum(exon_positions$length[c(1,2)]) + 1) > pos_tolerance), 
                            ((read_end_pos - sum(exon_positions$length[1]) + 1) > pos_tolerance))
    if(read_start_pass & read_end_pass)
    {
        write_to_sam_file(bam_data[[1]], i, output_sam)
    } else {
        read_out <- c(paste0("@", read_name), read_seq, "+", read_qual)
        fwrite(as.data.frame(read_out), output_fastq, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
}

#-- inputs --#
bam_file <- opt$input_bam
pos_tolerance <- opt$softclip
num_cores <- opt$threads
bam_chunk_size <- opt$chunk_size

bam_prefix <- tools::file_path_sans_ext(basename(bam_file))

exon_positions <- read.table(opt$exon_pos, header = FALSE, sep = "\t")
colnames(exon_positions) <- c("exon_id", "exon_start", "exon_end")
exon_positions <- as.data.table(exon_positions)
set(exon_positions, j = "exon_id", value = as.character(exon_positions$exon_id))
set(exon_positions, j = "exon_start", value = as.integer(exon_positions$exon_start))
set(exon_positions, j = "exon_end", value = as.integer(exon_positions$exon_end))
exon_positions$length <- exon_positions$exon_end - exon_positions$exon_start + 1

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

output_sam <- paste0(bam_prefix, ".filtered.sam")
output_fastq <- paste0(bam_prefix, ".wrongmap.fastq")
if(file.exists(output_sam)) invisible(file.remove(output_sam))
if(file.exists(output_fastq)) invisible(file.remove(output_fastq))

#-- processing --#
header <- scanBamHeader(bam_file)[[1]]$text
for(i in 1:length(header))
{
    if(names(header[i]) == "@HD" || names(header[i]) == "@SQ")
    {
        line <- c(names(header[i]), unlist(header[i]))
        fwrite(as.data.table(t(line)),
               output_sam,
               append = TRUE,
               quote = FALSE,
               sep = '\t',
               row.names = FALSE,
               col.names = FALSE)
    }
}

bam_fields <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual")
tag_fields <- c("NM", "AS", "XS")
param <- ScanBamParam(what = bam_fields, tag = tag_fields)

bam_read_handle <- open(BamFile(bam_file, yieldSize = bam_chunk_size))
bam_chunk_index <- 1
repeat
{
    bam_chunk <- scanBam(bam_read_handle, param = param, nThreads = num_cores)
    if(length(bam_chunk[[1]]$seq) == 0)
    {
        cat("done!\n")
        break
    } else {
        cat("\n")
        range_start <- format((bam_chunk_index - 1) * bam_chunk_size + 1, scientific = FALSE)
        range_end <- format((bam_chunk_index - 1) * bam_chunk_size + length(bam_chunk[[1]]$seq), scientific = FALSE)
        cat(paste0("processing chunk: ", range_start, " ~ ", range_end, "\n"))
        bam_chunk_index <- bam_chunk_index + 1
    }

    mclapply(1:length(bam_chunk[[1]]$seq), 
             function(i) process_read(i, bam_chunk), 
             mc.cores = num_cores)
}

close(bam_read_handle)
