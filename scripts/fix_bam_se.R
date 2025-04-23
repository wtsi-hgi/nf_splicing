#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Rsamtools", "Biostrings", "parallel", "optparse", "vroom")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bam"),          type = "character",  help = "input bam"),
                    make_option(c("-a", "--input_association"),  type = "character",  help = "barcode association file",  default = NULL),
                    make_option(c("-r", "--input_ref"),          type = "character",  help = "reference fasta file",      default = NULL),
                    make_option(c("-p", "--barcode_template"),   type = "character",  help = "barcode template",          default = "NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN"),
                    make_option(c("-m", "--barcode_marker"),     type = "character",  help = "barcode marker",            default = "CTACTGATTCGATGCAAGCTTG"),
                    make_option(c("-o", "--output_dir"),         type = "character",  help = "output directory",          default = getwd()),
                    make_option(c("-t", "--threads"),            type = "integer",    help = "number of threads",         default = 64),
                    make_option(c("-c", "--chunk_size"),         type = "integer",    help = "chunk size",                default = 10000),
                    make_option(c("-s", "--spliced"),            type = "logical",    help = "create spliced products",   default = FALSE,    action = "store_true"),
                    make_option(c("-l", "--library"),            type = "character",  help = "library type",              default = "muta"))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if(is.null(opt$input_bam)) stop("-b, input bam file is required!", call. = FALSE)
if(is.null(opt$input_association)) stop("-a, input association file is required!", call. = FALSE)
if(is.null(opt$input_ref)) stop("-r, input reference fasta file is required!", call. = FALSE)

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
                  paste0("XS:i:", bam_data$tag$XS[index]),
                  paste0("NH:i:", bam_data$tag$NH[index]),
                  paste0("MD:i:", bam_data$tag$MD[index]))

    fwrite(as.data.frame(t(read_out)), output_sam, append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

process_read <- function(i, bam_data, barcode_marker, barcode_template, barcode_length, barcode_map, ref_fasta, variant_map)
{
    # progress
    if(i %% 1000 == 0)
    {
        cat("\r", paste0("processed: ", i, "/", length(bam_data[[1]]$seq)))
    }

    # get the original read mapping info
    ref_id <- bam_data[[1]]$rname[i]
    ref_pos <- bam_data[[1]]$pos[i]
    cigar_lengths <- as.numeric(unlist(str_extract_all(bam_data[[1]]$cigar[i], "\\d+")))
    cigar_ops <- unlist(str_extract_all(bam_data[[1]]$cigar[i], "[A-Z]"))

    # fix the variant id by associated barcodes in the bam
    read_seq <- toupper(data.frame(bam_data[[1]]$seq[i])[1,1])
    barcode_seq <- ifelse(str_detect(barcode_template, "^[N]+$"),
                          get_barcode_by_marker(read_seq, barcode_marker, barcode_length),
                          get_barcode_by_template(read_seq, barcode_marker, barcode_template))

    if(!(barcode_seq %in% c("no barcode", "duplicate barcode")))
    {
        if(!is.na(barcode_map[barcode_seq])) 
        {
            bam_data[[1]]$rname[i] <- barcode_map[barcode_seq]
            output_line <- data.frame(barcode_map[barcode_seq], barcode_seq)
        } else {
            output_line <- data.frame("NA", barcode_seq)
        }
        fwrite(output_line, output_barcode, append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
    }

    write_to_sam_file(bam_data[[1]], i, output_sam)
    
    if(opt$spliced)
    {
        getseq_id <- ifelse(opt$library == "muta", ref_id, bam_data[[1]]$rname[i])

        spliced_seq <- DNAStringSet()
        for(j in seq_along(cigar_ops))
        {
            op <- cigar_ops[j]
            length <- cigar_lengths[j]

            if(op == "M" || op == "D")
            {
                ref_segment <- getSeq(ref_fasta, GRanges(seqnames = getseq_id, ranges = IRanges(ref_pos, ref_pos + length - 1)))
                spliced_seq <- c(spliced_seq, ref_segment)
                ref_pos <- ref_pos + length
            } else if(op == "N") {
                ref_pos <- ref_pos + length
            } else if (op == "S" || op == "I" || op == "H") {
                next
            }
        }

        output_id <- bam_data[[1]]$rname[i]
        output_line <- data.frame(output_id, variant_map[output_id], as.character(unlist(spliced_seq)))
        fwrite(output_line, output_tmp, append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
    }
}

#-- inputs --#
bam_file <- opt$input_bam
association_file <- opt$input_association
reference_file <- opt$input_ref
barcode_template <- opt$barcode_template
barcode_marker <- opt$barcode_marker
num_cores <- opt$threads
bam_chunk_size <- opt$chunk_size

barcode_length <- nchar(barcode_template)
bam_prefix <- tools::file_path_sans_ext(basename(opt$input_bam))

barcode_var <- fread(association_file, header = T, sep = '\t')
barcode_map <- setNames(barcode_var$varid, barcode_var$barcode)
variant_map <- setNames(barcode_var$variant, barcode_var$varid)
variant_map <- variant_map[!duplicated(variant_map)]

ref_fasta <- FaFile(reference_file)

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

output_barcode <- paste0(bam_prefix, ".barcodes.txt")
if(file.exists(output_barcode)) invisible(file.remove(output_file))

output_sam <- paste0(bam_prefix, ".fixed.sam")
if(file.exists(output_sam)) invisible(file.remove(output_sam))

if(opt$spliced)
{
    output_tmp <- paste0(bam_prefix, ".spliced_tmp.txt")
    if(file.exists(output_tmp)) invisible(file.remove(output_tmp))

    output_product <- paste0(bam_prefix, ".spliced_products.txt")
    if(file.exists(output_product)) invisible(file.remove(output_product))
}

#-- processing --#
header <- scanBamHeader(bam_file)[[1]]$text

ref_length <- 0
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

    if(names(header[i]) == "@SQ")
    {
        ref_length <- as.integer(unlist(strsplit(unlist(header[i])[2], ":"))[2])
    }
}

if(opt$library == "muta")
{
    variant_map_names <- names(variant_map)
    for(i in 1:length(variant_map_names))
    {
        line <- c("@SQ", paste0("SN:", variant_map_names[i]), paste0("LN:", ref_length))
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
open(ref_fasta)
bam_chunk_index <- 1
repeat
{
    bam_chunk <- scanBam(bam_read_handle, param = param, nThreads = num_cores)

    if(opt$library == "muta")
    {
        levels(bam_chunk[[1]]$rname) <- c(levels(bam_chunk[[1]]$rname), names(variant_map))
    }

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

    mclapply(1:length(bam_chunk[[1]]$seq), 
             function(i) process_read(i, bam_chunk, barcode_marker, barcode_template, barcode_length, barcode_map, ref_fasta, variant_map), 
             mc.cores = num_cores)
}

close(ref_fasta)
close(bam_read_handle)

if(opt$spliced)
{
    spliced_product <- as.data.table(vroom(output_tmp, delim = "\t", comment = "#", col_names = FALSE, show_col_types = FALSE))
    colnames(spliced_product) <- c("varid", "variant", "product")
    spliced_product_counts <- spliced_product[, .N, by = .(varid, variant, product)]
    colnames(spliced_product_counts)[3] <- "count"

    fwrite(spliced_product_counts, output_product, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    file.remove(output_tmp)
}
