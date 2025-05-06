#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Biostrings", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-b", "--barcode_association"), type = "character", help = "barcode association file",       default = NULL),
                    make_option(c("-c", "--canonical_barcodes"),  type = "character", help = "canonical barcode file",         default = NULL),
                    make_option(c("-n", "--novel_junctions"),     type = "character", help = "novel classified junction file", default = NULL),
                    make_option(c("-o", "--output_dir"),          type = "character", help = "output directory",               default = getwd()),
                    make_option(c("-p", "--prefix"),              type = "character", help = "output prefix",                  default = NULL))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
  print_help(opt_parser)
  quit(status = 1)
}

# Check if required arguments are provided
if(is.null(opt$barcode_association)) stop("-a, barcode association file is required!", call. = FALSE)
if(is.null(opt$canonical_barcodes))  stop("-b, canonical barcode file is required!", call. = FALSE)
if(is.null(opt$novel_junctions))  stop("-j, novel classified junction file is required!", call. = FALSE)

#-- inputs --#
barcode_association <- fread(opt$barcode_association, header = TRUE, sep = "\t")
var_idmap <- unique(barcode_association[, c("varid", "variant")])
setnames(var_idmap, c("varid", "variant"))

skipping_res <- fread(opt$canonical_barcodes, header = TRUE, sep = "\t")
novel_res <- fread(opt$novel_junctions, header = TRUE, sep = "\t")

output_prefix <- ifelse(is.null(opt$prefix), 
                        tools::file_path_sans_ext(basename(opt$novel_junctions)), 
                        opt$prefix)

categories <- c("intron_retention_I1",
                "intron_retention_I2",
                "intron_retension_5p_I1",
                "intron_retension_5p_I2",
                "intron_retension_3p_I1",
                "intron_retension_3p_I2",
                "exon_splicing_3p_E1",
                "exon_splicing_3p_E2",
                "exon_splicing_3p_E3",
                "exon_splicing_5p_E1",
                "exon_splicing_5p_E2",
                "exon_splicing_5p_E3",
                "exon_skipping_E2")

#-- outputs --#
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

#-- processing --#
# 1. initialize
splicing_matrix <- data.table(matrix(0L, nrow = dim(var_idmap)[1], ncol = 3))
setnames(splicing_matrix, c("varid", "canonical_inclusion_E2", "canonical_skipping_E2"))
splicing_matrix[, varid := as.character(varid)]
splicing_matrix$varid <- var_idmap$varid

# 2. integrate canonical events
skipping_res <- skipping_res[!is.na(varid)]
skipping_res_aggregated <- skipping_res[, .(count = sum(.SD$count)), by = .(varid, name)]

splicing_matrix_tmp <- merge(splicing_matrix, skipping_res_aggregated, by.x = "varid", by.y = "varid", all.x = TRUE)
splicing_matrix_tmp[, canonical_inclusion_E2 := ifelse(name == "E1_E2_E3", count, canonical_inclusion_E2)]
splicing_matrix_tmp[, canonical_skipping_E2 := ifelse(name == "E1_E3", count, canonical_skipping_E2)]
splicing_matrix_canonical <- splicing_matrix_tmp[, .(canonical_inclusion_E2 = sum(canonical_inclusion_E2, na.rm = TRUE),
                                                     canonical_skipping_E2 = sum(canonical_skipping_E2, na.rm = TRUE)), by = varid]

# 3. integrate novel events
novel_res <- novel_res[Annotation != "out_of_range"]
novel_res_split <- novel_res[, .(Annotation = unlist(strsplit(Annotation, ";")), 
                                 Cov = rep(Cov, lengths(strsplit(Annotation, ";")))), by = VarID]
novel_res_split <- novel_res_split[, .(Cov = sum(Cov)), by = .(VarID, Annotation)]
novel_res_wide <- dcast(novel_res_split, VarID ~ Annotation, value.var = "Cov", fun.aggregate = sum, fill = 0)

missing_categories <- setdiff(categories, names(novel_res_wide))
for(miss_cat in missing_categories)
{
    novel_res_wide[, (miss_cat) := as.integer(0)]
}

setcolorder(novel_res_wide, c("VarID", categories))

# 4. integrate all together
splicing_matrix_final <- merge(splicing_matrix_canonical, novel_res_wide, by.x = "varid", by.y = "VarID", all.x = TRUE)
splicing_matrix_final[is.na(splicing_matrix_final)] <- 0
setcolorder(splicing_matrix_final, c("varid", "canonical_inclusion_E2", "canonical_skipping_E2", categories))

fwrite(splicing_matrix_final, paste0(output_prefix, ".splicing_matrix.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
