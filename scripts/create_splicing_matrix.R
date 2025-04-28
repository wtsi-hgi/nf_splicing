#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("tidyverse", "data.table", "Biostrings", "optparse")
invisible(lapply(packages, quiet_library))

#-- options --#
option_list <- list(make_option(c("-a", "--input_association"),  type = "character", help = "barcode association file",       default = NULL),
                    make_option(c("-b", "--input_barcode"),      type = "character", help = "canonical barcode file",         default = NULL),
                    make_option(c("-j", "--input_junction"),     type = "character", help = "novel classified junction file", default = NULL),
                    make_option(c("-o", "--output_dir"),         type = "character", help = "output directory",               default = getwd()))
# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if(is.null(opt$input_bed)) stop("-b, input bed file is required!", call. = FALSE)
if(is.null(opt$exon_pos))  stop("-e, exon position file is required!", call. = FALSE)















library(tidyverse)
library(data.table)
library(Biostrings)

sample_id <- "HEK293T_Rep1"
associate_file <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/refseqs/varid_barcodes.unique.txt"
filtering_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/2_filtering"
mapping_dir <- "/lustre/scratch127/humgen/teams/hgi/fs18/splicing_analysis/benchmark/RON/3_mapping"

# reading files
skipping_file <- paste0(filtering_dir, "/", sample_id, "/", sample_id, ".barcodes_var.txt")
skipping_res <- as.data.table(read.table(skipping_file, header = T, sep = "\t"))

novel_file <- paste0(mapping_dir, "/", sample_id, "/", sample_id, ".junctions.txt")
novel_res <- as.data.table(read.table(novel_file, header = T, sep = "\t"))

barcode_association <- as.data.table(read.table(associate_file, header = TRUE, sep = "\t"))
var_idmap <- unique(barcode_association[, c("varid", "exon")])
setnames(var_idmap, c("varid", "vars"))

novel_res <- merge(novel_res, var_idmap, by.x = "VarID", by.y = "varid", all.x = TRUE)

# initialize
all_res <- data.table(matrix(0L, nrow = dim(var_idmap)[1], ncol = 3))
setnames(all_res, c("varid", "canonical_inclusion_E6", "canonical_skipping_E6"))
all_res[, varid := as.character(varid)]
all_res$varid <- var_idmap$varid

# integrate canonical events
skipping_res <- skipping_res[!is.na(varid)]
skipping_res_aggregated <- skipping_res[, .(count = sum(.SD$count)), by = .(varid, name)]

all_res_tmp <- merge(all_res, skipping_res_aggregated, by.x = "varid", by.y = "varid", all.x = TRUE)
all_res_tmp[, canonical_inclusion_E6 := ifelse(name == "e10_e11_e12", count, canonical_inclusion_E6)]
all_res_tmp[, canonical_skipping_E6 := ifelse(name == "e10_e12", count, canonical_skipping_E6)]
all_res_canonical <- all_res_tmp[, .(canonical_inclusion_E6 = sum(canonical_inclusion_E6, na.rm = TRUE),
                                     canonical_skipping_E6 = sum(canonical_skipping_E6, na.rm = TRUE)), by = varid]

# integrate novel events
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

# integrate all together
all_res_final <- merge(all_res_canonical, novel_res_wide, by.x = "varid", by.y = "VarID", all.x = TRUE)
all_res_final[is.na(all_res_final)] <- 0
setcolorder(all_res_final, c("varid", "canonical_inclusion_E6", "canonical_skipping_E6", categories))

write.table(all_res_final, paste0(sample_id, ".allres_by_var.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
