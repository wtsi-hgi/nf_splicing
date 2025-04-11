library(Rsamtools)
library(tidyverse)
library(data.table)
library(Biostrings)
library(parallel)

bam_dir <- "input_dir"
bam_file <- "input_bam"

barcode_template <- "NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN"
barcode_length <- nchar(barcode_template)
barcode_marker <- "CTACTGATTCGATGCAAGCTTG"

bam_prefix <- strsplit(bam_file, ".", fixed = T)[[1]][1]
setwd(bam_dir)

get_barcode_by_marker <- function(read_seq, barcode_marker, barcode_length)
{
    read_seq <- DNAString(read_seq)
    marker_match <- matchPattern(barcode_marker, read_seq, max.mismatch = 0, fixed = TRUE)

    if(length(marker_match) == 1) 
    {
        marker_start <- start(marker_match)[1]
        return(ifelse(marker_start > barcode_length, as.character(substr(read_seq, marker_start - barcode_length, marker_start - 1)), "no barcode"))
    } else if (length(marker_match) == 0) {
        marker_match_mis <- matchPattern(barcode_marker, read_seq, max.mismatch = 1, fixed = TRUE)

        if(length(marker_match_mis) == 1)
        {
            marker_start <- start(marker_match_mis)[1]
            return(ifelse(marker_start > barcode_length, as.character(substr(read_seq, marker_start - barcode_length, marker_start - 1)), "no barcode"))
        } else if(length(marker_match_mis) == 0) {
            return("no barcode")
        } else {
            return("duplicate barcode")
        }
    } else {
        return("duplicate barcode")
    }
}


get_barcode_by_template <- function(read_seq, barcode_marker, barcode_template)
{
    read_seq <- DNAString(read_seq)
    barcode_match <- matchPattern(barcode_template, read_seq, max.mismatch = 0, fixed = FALSE)

    if(length(barcode_match) == 1) 
    {
        return(as.character(barcode_match))
    } else if (length(barcode_match) == 0){
        # when barcode_match at the last part of read_seq, the length of matches barcode may less than barcode length
        # that is caused by the mismatch, and lead to the error for as.character()
        # eg:
        # seq <- "AATTCCGG"
        # marker <- "CGGA"
        # as.character(matchPattern(marker, seq, max.mismatch = 1, fixed = TRUE))        
        barcode_match_mis <- matchPattern(barcode_template, read_seq, max.mismatch = 1, fixed = FALSE, with.indels = TRUE)

        if(length(barcode_match_mis) == 1)
        {            
            return(as.character(barcode_match_mis))
        } else if(length(barcode_match_mis) == 0)
        {
            return("no barcode")
        } else {
            marker_match <- matchPattern(barcode_marker, read_seq, max.mismatch = 1, fixed = TRUE)
            if(length(marker_match) == 0)
            {
                return("no barcode")
            } else {
                marker_start <- marker_match@ranges@start[1]
                barcode_match_mis_ends <- end(barcode_match_mis@ranges)
                for(i in 1:length(barcode_match_mis_ends))
                {
                    if(marker_start == barcode_match_mis_ends[i] + 1)
                    {
                        return(as.character(barcode_match_mis)[i])
                    }
                }
                return("no barcode")
            }
        }
    } else {
        return("duplicate barcode")
    }
}

create_barcode_list <- function(i, bam_data, barcode_marker, barcode_template, barcode_length)
{
    if(i %% 5000 == 0)
    {
        cat("\r", paste0("processed: ", i, "/", length(bam_data[[1]]$seq)))
    }

    read_ref <- as.character(bam_data[[1]]$rname[i])
    if(!(read_ref %in% names(barcodes_chunk)))
    {
        barcodes_chunk[[read_ref]] <- vector()
    }

    read_seq <- toupper(data.frame(bam_data[[1]]$seq[i])[1,1])
    barcode_seq <- ifelse(str_detect(barcode_template, "^[N]+$"),
                          get_barcode_by_marker(read_seq, barcode_marker, barcode_length),
                          get_barcode_by_template(read_seq, barcode_marker, barcode_template))

    barcodes_chunk[[read_ref]] <- c(barcodes_chunk[[read_ref]], barcode_seq)

    return(barcodes_chunk)
}

# limit the reading info to save RAM, and also limit the yeild
param <- ScanBamParam(what = c("rname", "seq"))
#param <- ScanBamParam(what = scanBamWhat(), tag = c("NM", "AS", "XS"))

num_cores <- 32
bam_chunk_size <- 100000

# read bam into chunks
barcodes_list <- list()
bam_file_handle <- open(BamFile(bam_file, yieldSize = bam_chunk_size))
bam_chunk_index <- 1
repeat
{
    barcodes_chunk <- list()
    bam_chunk <- scanBam(bam_file_handle, param = param, nThreads = num_cores)
    if(length(bam_chunk[[1]]$seq) == 0)
    {
        cat("done!\n")
        break
    } else {
        cat("\n")
        range_start <- format((bam_chunk_index - 1) * bam_chunk_size + 1, scientific = FALSE)
        range_end <- format((bam_chunk_index - 1) * bam_chunk_size + length(bam_chunk[[1]]$seq), scientific = FALSE)
        cat(paste0("chunk range: ", range_start, " ~ ", range_end))
        bam_chunk_index <- bam_chunk_index + 1
    }

    barcodes_chunk <- mclapply(1:length(bam_chunk[[1]]$seq), 
                               function(i) create_barcode_list(i, bam_chunk, barcode_marker, barcode_template, barcode_length), 
                               mc.cores = num_cores)
    barcodes_chunk_list <- split(unlist(barcodes_chunk), names(unlist(barcodes_chunk)))

    for(name in names(barcodes_chunk_list))
    {
        if(name %in% names(barcodes_list))
        {
            barcodes_list[[name]] <- c(barcodes_list[[name]] , barcodes_chunk_list[[name]])
        } else {
            barcodes_list[[name]] <- barcodes_chunk_list[[name]]
        }
    }
}
close(bam_file_handle)

#bam_data <- scanBam(bam_file, param = param, nThreads = num_cores)
#barcodes <- list()
#barcodes <- mclapply(1:length(bam_data[[1]]$seq), function(i) create_barcode_list(i), mc.cores = num_cores)

#barcodes_list <- split(unlist(barcodes), names(unlist(barcodes)))
for(i in 1:length(barcodes_list))
{
    barcodes_list[[i]] <- as.data.table(table(barcodes_list[[i]]))
    colnames(barcodes_list[[i]]) <- c("barcode", "counts")

    barcodes_list[[i]]$name <- names(barcodes_list)[i]
    barcodes_list[[i]] <- barcodes_list[[i]][, c("name", "barcode", "counts")]
}

barcodes_out <- do.call(rbind, barcodes_list)
output_file <- paste0(bam_prefix, ".barcodes.txt")
write.table(barcodes_out, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
