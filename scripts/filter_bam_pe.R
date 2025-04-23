#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Rsamtools", "Biostrings", "parallel", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bam"),  type = "character", help = "input bam"),
                    make_option(c("-e", "--exon_pos"),   type = "character", help = "exon position file"),
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

revcomp <- function(seq) {
    seq <- toupper(seq)
    splits <- strsplit(seq, "")[[1]]
    reversed <- rev(splits)
    seq_rev <- paste(reversed, collapse = "")
    seq_rev_comp <- chartr("ATCG", "TAGC", seq_rev)
    return(seq_rev_comp)
}

calc_softclip_lens <- function(cigar)
{
    first_softclip <- 0
    last_softclip <- 0
    
    if(grepl("^[0-9]+S", cigar)) first_softclip <- as.numeric(sub("S", "", regmatches(cigar, regexpr("^[0-9]+S", cigar))))
    if(grepl("[0-9]+S$", cigar)) last_softclip <- as.numeric(sub("S", "", regmatches(cigar, regexpr("[0-9]+S$", cigar))))

    return(c(first_softclip, last_softclip))
}

process_read <- function(i, bam_data)
{
    read1_index <- i
    read2_index <- i + 1

    # R1 is Forward Strand: 99
    # R2 is Reverse Strand: 147
    # R1 is Reverse Strand: 83
    # R2 is Forward Strand: 163
    # for analysis, left (forward) is read1, right (reverse) is read2
    if (bam_data[[1]]$flag[read1_index] == 147 | bam_data[[1]]$flag[read1_index] == 83)
    {
        read1_index <- i + 1
        read2_index <- i
    }

    read1_name <- bam_data[[1]]$qname[read1_index]
    read1_start_pos <- bam_data[[1]]$pos[read1_index]
    read1_cigar <- bam_data[[1]]$cigar[read1_index]
    read1_end_pos <- calc_end_pos(read1_start_pos, read1_cigar)
    read1_ref <- bam_data[[1]]$rname[read1_index]
    read1_seq <- data.frame(bam_data[[1]]$seq[read1_index])[1,1]
    read1_qual <- data.frame(bam_data[[1]]$qual[read1_index])[1,1]

    read2_name <- bam_data[[1]]$qname[read2_index]
    read2_start_pos <- bam_data[[1]]$pos[read2_index]
    read2_cigar <- bam_data[[1]]$cigar[read2_index]
    read2_end_pos <- calc_end_pos(read2_start_pos, read2_cigar)
    read2_ref <- bam_data[[1]]$rname[read2_index]
    read2_seq <- revcomp(data.frame(bam_data[[1]]$seq[read2_index])[1,1])
    read2_qual <- data.frame(bam_data[[1]]$qual[read2_index])[1,1]

    if((read1_name == read2_name) & (read1_ref == read2_ref))
    {
        # read1 alignment should start at least 5 bp left-away from target exon start postion due to exon 5
        # read2 alignment should end near the end of reference due to barcode softclip
        # if the library is long, the read may not cover the end part of last exon, allow 5 bp difference
        read_start_pass <- ((exon_positions$exon_end[1] - read1_start_pos) > pos_tolerance)
        read_end_pass <- ifelse(read1_ref == "E1_E2_E3", 
                                ((read2_end_pos - sum(exon_positions$length[c(1,2)]) + 1) > pos_tolerance), 
                                ((read2_end_pos - sum(exon_positions$length[1]) + 1) > pos_tolerance))
        if(read_start_pass & read_end_pass)
        {
            # read1 should not have end softclip, if has, the read1 should end near the end of reference
            condition1 <- ifelse(str_detect(read1_cigar, "\\d+S$"),
                                 ifelse((calc_softclip_lens(read1_cigar)[2] > pos_tolerance), 
                                        ((chrs[read1_ref] - read1_end_pos) < pos_tolerance), 
                                        TRUE), 
                                 TRUE)

            # read2 should not have start softclip, if has, the length of softclip should be less than 5bp                     
            condition2 <- ifelse(str_detect(read2_cigar, "^\\d+S"),
                                 ifelse((calc_softclip_lens(read2_cigar)[1] <= pos_tolerance), 
                                        (read2_start_pos >= read1_end_pos), 
                                        FALSE), 
                                 (read2_start_pos >= read1_end_pos))

            if(condition1 & condition2)
            {
                write_to_sam_file(bam_data[[1]], read1_index, output_sam)
                write_to_sam_file(bam_data[[1]], read2_index, output_sam)
            } else {
                read1_out <- c(paste0("@", read1_name), read1_seq, "+", read1_qual)
                read2_out <- c(paste0("@", read2_name), read2_seq, "+", read2_qual)  

                fwrite(as.data.frame(read1_out), output_fastq_r1, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)     
                fwrite(as.data.frame(read2_out), output_fastq_r2, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) 
            }
        } else {
            read1_out <- c(paste0("@", read1_name), read1_seq, "+", read1_qual)
            read2_out <- c(paste0("@", read2_name), read2_seq, "+", read2_qual)  

            fwrite(as.data.frame(read1_out), output_fastq_r1, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)     
            fwrite(as.data.frame(read2_out), output_fastq_r2, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) 
        }
    } else {
        read1_out <- c(paste0("@", read1_name), read1_seq, "+", read1_qual)
        read2_out <- c(paste0("@", read2_name), read2_seq, "+", read2_qual)  

        fwrite(as.data.frame(read1_out), output_fastq_r1, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)     
        fwrite(as.data.frame(read2_out), output_fastq_r2, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) 
    }
}

#-- inputs --#
bam_file <- opt$input_bam
pos_tolerance <- opt$softclip
num_cores <- opt$threads
bam_chunk_size <- opt$chunk_size

bam_prefix <- tools::file_path_sans_ext(basename(opt$input_bam))

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
# need to check if wrongmap has the same read strand with unmap
output_fastq_r1 <- paste0(bam_prefix, ".wrongmap.r2.fastq")
output_fastq_r2 <- paste0(bam_prefix, ".wrongmap.r1.fastq")
if(file.exists(output_sam)) invisible(file.remove(output_sam))
if(file.exists(output_fastq_r1)) invisible(file.remove(output_fastq_r1))
if(file.exists(output_fastq_r2)) invisible(file.remove(output_fastq_r2))

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

    mclapply(seq(1,length(bam_chunk[[1]]$seq),2), 
             function(i) process_read(i, bam_chunk), 
             mc.cores = num_cores)
}

close(bam_read_handle)
