#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Rsamtools", "Biostrings", "parallel", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bam"),         type = "character", help = "input bam",                default = NULL),
                    make_option(c("-a", "--input_association"), type = "character", help = "barcode association file", default = NULL),
                    make_option(c("-p", "--barcode_template"),  type = "character", help = "barcode template",         default = "NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN"),
                    make_option(c("-m", "--barcode_marker"),    type = "character", help = "barcode marker",           default = "CTACTGATTCGATGCAAGCTTG"),
                    make_option(c("-o", "--output_dir"),        type = "character", help = "output directory",         default = getwd()),
                    make_option(c("-t", "--threads"),           type = "integer",   help = "number of threads",        default = 64),
                    make_option(c("-c", "--chunk_size"),        type = "integer",   help = "chunk size",               default = 10000))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
  print_help(opt_parser)
  quit(status = 1)
}

# Check if required arguments are provided
if(is.null(opt$input_bam)) stop("-b, input bam file is required!", call. = FALSE)

#-- function --#
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

# create_barcode_list <- function(i, bam_data, barcode_marker, barcode_template, barcode_length)
create_barcode_row <- function(i, bam_data, barcode_marker, barcode_template, barcode_length)
{
    if(i %% 5000 == 0)
    {
        cat("\r", paste0("processed: ", i, "/", length(bam_data[[1]]$seq)))
    }

    # different experiment needs to look at different flag
    # if reference is plus strand, look at 147
    # if reference is minus strand, look at 83    
    if(bam_data[[1]]$flag[i] == 83)
    {
        read_ref <- as.character(bam_data[[1]]$rname[i])
        read_seq <- toupper(data.frame(bam_data[[1]]$seq[i])[1,1])
        barcode_seq <- ifelse(str_detect(barcode_template, "^[N]+$"),
                              get_barcode_by_marker(read_seq, barcode_marker, barcode_length),
                              get_barcode_by_template(read_seq, barcode_marker, barcode_template))

        # if(!(read_ref %in% names(barcodes_chunk)))
        # {
        #     barcodes_chunk[[read_ref]] <- vector()
        # }
        # barcodes_chunk[[read_ref]] <- c(barcodes_chunk[[read_ref]], barcode_seq)
        return(data.table(ref_id = read_ref, barcode = barcode_seq))
    }

    # return(barcodes_chunk)
    return(data.table(ref_id = character(0), barcode = character(0)))
}

#-- inputs --#
bam_file <- opt$input_bam
num_cores <- opt$threads
bam_chunk_size <- opt$chunk_size
barcode_template <- opt$barcode_template
barcode_marker <- opt$barcode_marker

barcode_length <- nchar(barcode_template)
bam_prefix <- tools::file_path_sans_ext(basename(bam_file))

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

output_file <- paste0(bam_prefix, ".barcodes.txt")
if(file.exists(output_file)) invisible(file.remove(output_file))

#-- processing --#
param <- ScanBamParam(what = c("rname", "flag", "seq"))

barcodes_list <- list()
bam_file_handle <- open(BamFile(bam_file, yieldSize = bam_chunk_size))
bam_chunk_index <- 1
repeat
{
    # barcodes_chunk <- list()
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

    barcodes_chunk_rows <- mclapply(1:length(bam_chunk[[1]]$seq),
                                    function(i) create_barcode_row(i, bam_chunk, barcode_marker, barcode_template, barcode_length),
                                    mc.cores = num_cores)
    barcodes_chunk_dt <- rbindlist(barcodes_chunk_rows)
    barcodes_list[[length(barcodes_list) + 1]] <- barcodes_chunk_dt

    # barcodes_chunk <- mclapply(1:length(bam_chunk[[1]]$seq), 
    #                            function(i) create_barcode_list(i, bam_chunk, barcode_marker, barcode_template, barcode_length), 
    #                            mc.cores = num_cores)
    # barcodes_chunk_list <- split(unlist(barcodes_chunk), names(unlist(barcodes_chunk)))

    # for(name in names(barcodes_chunk_list))
    # {
    #     if(name %in% names(barcodes_list))
    #     {
    #         barcodes_list[[name]] <- c(barcodes_list[[name]] , barcodes_chunk_list[[name]])
    #     } else {
    #         barcodes_list[[name]] <- barcodes_chunk_list[[name]]
    #     }
    # }
}
close(bam_file_handle)

if(length(barcodes_list) == 0)
{
    fwrite(data.table(), output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else {
    # for(i in 1:length(barcodes_list))
    # {
    #     barcodes_list[[i]] <- as.data.table(table(barcodes_list[[i]]))
    #     colnames(barcodes_list[[i]]) <- c("barcode", "count")

    #     barcodes_list[[i]]$name <- names(barcodes_list)[i]
    #     barcodes_list[[i]] <- barcodes_list[[i]][, c("name", "barcode", "count")]
    # }

    # barcodes_out <- as.data.table(do.call(rbind, barcodes_list))

    barcodes_out <- rbindlist(barcodes_list)
    barcodes_out <- barcodes_out[, .N, by = .(ref_id, barcode)]
    setnames(barcodes_out, "N", "count")
    if(is.null(opt$input_association))
    {
        fwrite(barcodes_out, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    } else {
        barcode_variant <- fread(opt$input_association, sep = "\t", header = TRUE)
        barcodes_out <- merge(barcodes_out, barcode_variant[, .(barcode, varid)], by = "barcode", all.x = TRUE)
        barcodes_out$varid[is.na(barcodes_out$varid)] <- "NA"
        barcodes_out <- barcodes_out[, .(ref_id, barcode, varid, count)]
        setorder(barcodes_out, ref_id)
        fwrite(barcodes_out, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
}