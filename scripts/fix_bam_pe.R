#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Rsamtools", "Biostrings", "parallel", "optparse", "vroom")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--input_bam"),         type = "character", help = "input bam",                default = NULL),
                    make_option(c("-l", "--lib_type"),          type = "character", help = "library type",             default = "muta_exon"),
                    make_option(c("-a", "--input_association"), type = "character", help = "barcode association file", default = NULL),
                    make_option(c("-r", "--input_ref"),         type = "character", help = "reference fasta file",     default = NULL),
                    make_option(c("-p", "--barcode_template"),  type = "character", help = "barcode template",         default = "NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN"),
                    make_option(c("-m", "--barcode_marker"),    type = "character", help = "barcode marker",           default = "CTACTGATTCGATGCAAGCTTG"),
                    make_option(c("-o", "--output_dir"),        type = "character", help = "output directory",         default = getwd()),
                    make_option(c("-t", "--threads"),           type = "integer",   help = "number of threads",        default = 64),
                    make_option(c("-c", "--chunk_size"),        type = "integer",   help = "chunk size",               default = 10000),
                    make_option(c("-s", "--spliced"),           type = "logical",   help = "create spliced products",  default = FALSE,    action = "store_true"))
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
if(is.null(opt$lib_type))  stop("-l, library type is required!", call. = FALSE)
if(is.null(opt$input_association)) stop("-a, input association file is required!", call. = FALSE)
if(is.null(opt$input_ref)) stop("-r, input reference fasta file is required!", call. = FALSE)

valid_lib_types <- c("random_intron", "random_exon", "muta_exon")
if(!(opt$lib_type %in% valid_lib_types))
{
    stop(paste0("-l, library type must be one of: ", paste(valid_lib_types, collapse = ", ")), call. = FALSE)
}

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
                  ifelse(is.na(bam_data$mrnm[index]), "=", as.character(bam_data$mrnm[index])),
                  ifelse(is.na(bam_data$mpos[index]), 0, bam_data$mpos[index]),
                  bam_data$isize[index],
                  data.frame(bam_data$seq[index])[1,1],
                  data.frame(bam_data$qual[index])[1,1],
                  paste0("NM:i:", bam_data$tag$NM[index]),
                  paste0("AS:i:", bam_data$tag$AS[index]),
                  paste0("NH:i:", bam_data$tag$NH[index]),
                  paste0("MD:i:", bam_data$tag$MD[index]))

    fwrite(as.data.frame(t(read_out)), output_sam, append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}


process_read <- function(i, bam_data, barcode_marker, barcode_template, barcode_length, barcode_map, ref_fasta, variant_map)
{
    # R1 is Forward Strand: 99
    # R2 is Reverse Strand: 147
    # R1 is Reverse Strand: 83
    # R2 is Forward Strand: 163
    # for analysis, left (forward) is read1, right (reverse) is read2 with barcode
    if(bam_data[[1]]$flag[i] == 147 | bam_data[[1]]$flag[i] == 83)
    {
        read1_index <- i + 1
        read2_index <- i
    } else {
        read1_index <- i
        read2_index <- i + 1
    }

    # get the original read mapping info
    read1_ref <- bam_data[[1]]$rname[read1_index]
    read2_ref <- bam_data[[1]]$rname[read2_index]
    read1_pos <- bam_data[[1]]$pos[read1_index]
    read2_pos <- bam_data[[1]]$pos[read2_index]
    read1_cigar_ops <- unlist(str_extract_all(bam_data[[1]]$cigar[read1_index], "[A-Z]"))
    read2_cigar_ops <- unlist(str_extract_all(bam_data[[1]]$cigar[read2_index], "[A-Z]"))
    read1_cigar_lengths <- as.numeric(unlist(str_extract_all(bam_data[[1]]$cigar[read1_index], "\\d+")))
    read2_cigar_lengths <- as.numeric(unlist(str_extract_all(bam_data[[1]]$cigar[read2_index], "\\d+")))
    read1_seq <- data.frame(bam_data[[1]]$seq[read1_index])[1,1]
    read2_seq <- data.frame(bam_data[[1]]$seq[read2_index])[1,1]
    
    # fix the variant id by associated barcodes in the bam
    read_seq <- toupper(read2_seq)
    barcode_seq <- ifelse(str_detect(barcode_template, "^[N]+$"),
                          get_barcode_by_marker(read_seq, barcode_marker, barcode_length),
                          get_barcode_by_template(read_seq, barcode_marker, barcode_template))

    if(!is.na(barcode_map[barcode_seq])) 
    {
        bam_data[[1]]$rname[read1_index] <- barcode_map[barcode_seq]
        bam_data[[1]]$rname[read2_index] <- barcode_map[barcode_seq]
        bam_data[[1]]$mrnm[read1_index] <- barcode_map[barcode_seq]
        bam_data[[1]]$mrnm[read2_index] <- barcode_map[barcode_seq]
        dt_out <- data.table(barcode = barcode_seq, varid = barcode_map[barcode_seq])
    } else {
        dt_out <- data.table(barcode = barcode_seq, varid = "NA")
    }

    write_to_sam_file(bam_data[[1]], read1_index, output_sam)
    write_to_sam_file(bam_data[[1]], read2_index, output_sam)

    if(opt$spliced)
    {
        # random exon/intron library, its hisat2 ref includes all the random intron sequences, barcode per random
        # muta exon library, its hisat2 ref includes all the wild-type sequences, barcode per variant
        getseq_id <- ifelse(opt$lib_type == "muta_exon", as.character(read1_ref), as.character(bam_data[[1]]$rname[read1_index]))
        
        read1_spliced_seq <- DNAStringSet()

        for(j in seq_along(read1_cigar_ops))
        {
            op <- read1_cigar_ops[j]
            length <- read1_cigar_lengths[j]

            if(op == "M" || op == "D")
            {
                ref_segment <- subseq(ref_sequences[[getseq_id]], start = ref_pos, end = ref_pos + length - 1)
                spliced_seq <- c(spliced_seq, DNAStringSet(ref_segment))
                read1_pos <- read1_pos + length
            } else if(op == "N") {
                read1_pos <- read1_pos + length
            } else if (op == "S" || op == "I" || op == "H") {
                next
            }
        }

        read2_spliced_seq <- DNAStringSet()
        for(j in seq_along(read2_cigar_ops))
        {
            op <- read2_cigar_ops[j]
            length <- read2_cigar_lengths[j]

            if(op == "M" || op == "D")
            {
                ref_segment <- getSeq(ref_fasta, GRanges(seqnames = getseq_id, ranges = IRanges(read2_pos, read2_pos + length - 1)))
                read2_spliced_seq <- c(read2_spliced_seq, ref_segment)
                read2_pos <- read2_pos + length
            } else if(op == "N") {
                read2_pos <- read2_pos + length
            } else if (op == "S" || op == "I" || op == "H") {
                next
            }
        }

        # read1_pos and read2_pos have been updated to the end of alignments, so use bam_data
        # get the gap sequence between read1 and read2
        # normally read1 and read2 should not overlap in the middle as reads have been separated
        if(read1_pos < bam_data[[1]]$pos[read2_index])
        {
            read_gap <- getSeq(ref_fasta, GRanges(seqnames = getseq_id, ranges = IRanges(read1_pos, bam_data[[1]]$pos[read2_index] - 1)))
            spliced_seq <- paste0(as.character(unlist(read1_spliced_seq)), as.character(read_gap), as.character(unlist(read2_spliced_seq)))
        } else {
            spliced_seq <- paste0(as.character(unlist(read1_spliced_seq)), as.character(unlist(read2_spliced_seq)))
        }

        output_id <- bam_data[[1]]$rname[i]
        output_line <- ifelse(opt$lib_type == "muta_exon", 
                              data.frame(output_id, getseq_id, spliced_seq)
                              data.frame(output_id, variant_map[output_id], spliced_seq))
        # output_line <- data.frame(output_id, variant_map[output_id], spliced_seq)
        fwrite(output_line, output_tmp, append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
    }

    return(dt_out)
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
bam_prefix <- tools::file_path_sans_ext(basename(bam_file))

barcode_var <- fread(association_file, header = T, sep = '\t')
barcode_map <- setNames(barcode_var$varid, barcode_var$barcode)
variant_map <- setNames(barcode_var$variant, barcode_var$varid)
variant_map <- variant_map[!duplicated(variant_map)]

ref_sequences <- readDNAStringSet(reference_file)

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

output_barcode <- paste0(bam_prefix, ".barcodes.txt")
if(file.exists(output_barcode)) invisible(file.remove(output_barcode))

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

if(opt$lib_type == "muta_exon")
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
tag_fields <- c("NM", "AS", "NH", "MD")
param <- ScanBamParam(what = bam_fields, tag = tag_fields)

dt_barcodes <- data.table(barcode = character(), varid = character()) 
bam_read_handle <- open(BamFile(bam_file, yieldSize = bam_chunk_size))
bam_chunk_index <- 1
repeat
{
    chunk_results <- list()
    bam_chunk <- scanBam(bam_read_handle, param = param, nThreads = num_cores)

    if(opt$lib_type == "muta_exon")
    {
        levels(bam_chunk[[1]]$rname) <- c(levels(bam_chunk[[1]]$rname), names(variant_map))
        levels(bam_chunk[[1]]$mrnm) <- c(levels(bam_chunk[[1]]$mrnm), names(variant_map))
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

    # mclapply cannot modify global variables, so we need to pass all the variables to the function
    # and then collect the results in a list
    chunk_results <- mclapply(seq(1,length(bam_chunk[[1]]$seq),2),
                              function(i) process_read(i, bam_chunk, barcode_marker, barcode_template, barcode_length, barcode_map, ref_fasta, variant_map), 
                              mc.cores = num_cores)
    dt_barcodes <- rbind(dt_barcodes, rbindlist(chunk_results))
}

close(bam_read_handle)

dt_barcodes <- dt_barcodes[, .N, by = .(barcode, varid)]
colnames(dt_barcodes) <- c("barcode", "varid", "count")
fwrite(dt_barcodes, output_barcode, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

if(opt$spliced)
{
    spliced_product <- as.data.table(vroom(output_tmp, delim = "\t", comment = "#", col_names = FALSE, show_col_types = FALSE))
    colnames(spliced_product) <- c("varid", "variant", "product")
    spliced_product_counts <- spliced_product[, .N, by = .(varid, variant, product)]
    colnames(spliced_product_counts)[3] <- "count"

    fwrite(spliced_product_counts, output_product, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    file.remove(output_tmp)
}
