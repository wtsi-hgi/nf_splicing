create_barplots <- function(barcode_association, sample_reps, total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads)
{
    summary_reads <- as.data.table(cbind(sample_reps, total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads))
    colnames(summary_reads) <- c("sample_reps", "total_reads", "merged_reads", "unmerged_reads", "inclusion_reads", "skipping_reads", "map_reads", "unexplain_reads")

    summary_reads <- summary_reads %>% mutate(across(c(total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads), as.numeric))
    summary_pct <- summary_reads %>%
                    mutate(pct_merged = round((merged_reads / total_reads) * 100, 1),
                           pct_unmerged = round((unmerged_reads / total_reads) * 100, 1),
                           library_coverage = as.integer(total_reads / length(unique(barcode_association$var_id))),
                           pct_inclusion = round((inclusion_reads / total_reads) * 100, 1),
                           pct_skipping = round((skipping_reads / total_reads) * 100, 1),
                           pct_map = round((map_reads / total_reads) * 100, 1),
                           pct_unexplain = round((unexplain_reads / total_reads) * 100, 1))

    tib_pct_total <- summary_pct %>% 
                        select(sample_reps, pct_merged, pct_unmerged) %>%
                        pivot_longer(cols = -sample_reps, names_to = "type", values_to = "pct")
    tib_pct_total$type <- factor(tib_pct_total$type, levels = c("pct_unmerged", "pct_merged"))
    tib_cov <- summary_pct %>% select(sample_reps, library_coverage)
    tib_cov$type <- "coverage"

    tib_pct_canonical <- summary_pct %>% 
                            select(sample_reps, pct_inclusion, pct_skipping) %>%
                            pivot_longer(cols = -sample_reps, names_to = "type", values_to = "pct")

    tib_pct_novel <- summary_pct %>% 
                        select(sample_reps, pct_map, pct_unexplain) %>%
                        pivot_longer(cols = -sample_reps, names_to = "type", values_to = "pct")

    select_colors <- select_colorblind("col15")[1:2]
    fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

    p1 <- ggplot(tib_pct_total,  aes(x = sample_reps, y = pct, fill = type)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_fill_manual(values = fill_colors) +
            labs(x = NULL, y = "percent") +
            theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank()) +
            theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
            theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
            theme(axis.text = element_text(size = 8, face = "bold")) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)

    p2 <- ggplot(tib_pct_canonical,  aes(x = sample_reps, y = pct, fill = type)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_fill_manual(values = fill_colors) +
            labs(x = NULL, y = "percent") +
            theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank()) +
            theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
            theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
            theme(axis.text = element_text(size = 8, face = "bold")) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)

    p3 <- ggplot(tib_pct_novel,  aes(x = sample_reps, y = pct, fill = type)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_fill_manual(values = fill_colors) +
            labs(x = NULL, y = "percent") +
            theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank()) +
            theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
            theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
            theme(axis.text = element_text(size = 8, face = "bold")) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)
    
    return(list(p1, p2, p3))
}

create_venn_diagrams <- function(canonical_barcodes, novel_barcodes, barcode_association)
{
    list_plots <- list()
    for(i in seq_along(canonical_barcodes))
    {
        venn_list <- list(unique(canonical_barcodes[[i]]$barcode), unique(novel_barcodes[[i]]$barcode), barcode_association$barcode)
        names(venn_list) <- c("canonical", "novel", "association")
        p <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.2) + scale_fill_gradient(low = "ivory", high = "tomato")
        list_plots[[i]] <- p
    }
    return(list_plots)
}

create_junction_plots <- function(classified_junctions)
{
    list_junctions <- lapply(classified_junctions, function(dt) { paste(dt$chrom, dt$donor, dt$acceptor, sep = "_") })
    p <- ggVennDiagram(list_junctions, label_alpha = 0, edge_size = 0.2) + scale_fill_gradient(low = "ivory", high = "tomato")

    dt_list <- lapply(seq_along(classified_junctions), function(i) {
                    dt <- classified_junctions[[i]][, .(chrom, donor, acceptor, annotation, coverage)]
                    setnames(dt, "coverage", sample_reps[i])
                    dt })
    dt_junctions <- Reduce(function(x, y) merge(x, y, by = c("chrom", "donor", "acceptor", "annotation"), all = TRUE), dt_list)
    
    cov_cols <- setdiff(names(dt_junctions), c("chrom", "donor", "acceptor", "annotation"))
    dt_junctions <- dt_junctions[rowSums(!is.na(dt_junctions[, ..cov_cols])) >= 2]
    dt_junctions[, avg_cov := rowMeans(.SD, na.rm = TRUE), .SDcols = cov_cols]
    setnames(dt_junctions, "chrom", "var_id")

    return(list(p, dt_junctions))
}

rescale_junctions <- function(junctions, exons, lib_type)
{   
    exon_template <- exons[1:3,]
    introns_template <- exon_template[, .(intron_id   = paste0("I", .I - 1),
                                          intron_start = shift(exon_end, type = "lag") + 1,
                                          intron_end   = exon_start - 1), by = .(var_id)][!is.na(intron_start)]

    if(lib_type == "random_intron" || lib_type == "random_exon")
    {
        return(list(junctions, exon_template, introns_template))
    }
    
    junctions_tmp <- copy(junctions)
    junctions_tmp[, var_id := sub("(/.*)$", "", var_id)]

    if(lib_type == "muta_exon")
    {
        targets <- copy(exons)
        setnames(targets, c("exon_id", "exon_start", "exon_end"), c("target_id", "target_start", "target_end"))
        targets <- merge(targets, 
                         exon_template[, .(target_id = exon_id, template_start = exon_start, template_end = exon_end)], 
                         by = "target_id", 
                         all.x = TRUE)

    } else {
        targets <- exons[, .(target_id    = paste0("I", seq_len(.N) - 1),
                             target_start = shift(exon_end, type = "lag") + 1,
                             target_end   = exon_start - 1), by = var_id][!is.na(target_start)]
        targets <- merge(targets, 
                         introns_template[, .(target_id = intron_id, template_start = intron_start, template_end = intron_end)], 
                         by = "target_id", 
                         all.x = TRUE)
    }
    
    junctions_tmp <- merge(junctions_tmp, targets, by = "var_id", allow.cartesian = TRUE)

    junctions_tmp[, donor_rescale := fifelse(donor >= target_start & donor <= target_end,
                             round(template_start + (donor - target_start) / (target_end - target_start) * (template_end - template_start)),
                             donor)]

    junctions_tmp[, acceptor_rescale := fifelse(acceptor >= target_start & acceptor <= target_end,
                             round(template_start + (acceptor - target_start) / (target_end - target_start) * (template_end - template_start)),
                             acceptor)]

    junctions_tmp <- junctions_tmp[, .(donor_rescale = min(donor_rescale), acceptor_rescale = min(acceptor_rescale)), 
                                    by = .(var_id, donor, acceptor, annotation, I4_Rep1, I4_Rep2, I4_Rep3, avg_cov)]
    junctions_tmp <- junctions_tmp %>% select(-c(donor, acceptor))
    setnames(junctions_tmp, c("donor_rescale", "acceptor_rescale"), c("donor", "acceptor"))
    setcolorder(junctions_tmp, colnames(junctions))

    return(list(junctions, exon_template, introns_template))
}

create_junction_distribution <- function(junctions, exons, introns)
{
    junctions_a1b10 <- junctions %>% filter(avg_cov >= 1 & avg_cov < 10)
    junctions_a10b100 <- junctions %>% filter(avg_cov >= 10 & avg_cov < 100)
    junctions_a100b1000 <- junctions %>% filter(avg_cov >= 100 & avg_cov < 1000)
    junctions_a1000b10000 <- junctions %>% filter(avg_cov >= 1000 & avg_cov < 10000)
    junctions_a10000 <- junctions %>% filter(avg_cov >= 10000)

    list_junctions <- list(junctions, junctions_a1b10, junctions_a10b100, junctions_a100b1000, junctions_a1000b10000, junctions_a10000)
    names(list_junctions) <- c("avg_cov:all", "avg_cov:1-10", "avg_cov:10-100", "avg_cov:100-1000", "avg_cov:1000-10000", "avg_cov:10000+")
    list_diagramplots <- list()
    list_scatterplots <- list()

    for(i in seq_along(list_junctions))
    {
        p1 <- ggplot() +
                geom_rect(data = exons, fill = t_col("darkgreen", 0.6), color = t_col("darkgreen", 0.6), aes(xmin = exon_start, xmax = exon_end, ymin = 0.8, ymax = 1)) +
                geom_rect(data = introns, fill = "grey", color = "grey", aes(xmin = intron_start, xmax = intron_end, ymin = 0.87, ymax = 0.93)) +
                geom_curve(data = list_junctions[[i]], aes(x = donor, xend = acceptor, y = 1, yend = 1, color = avg_cov, alpha = avg_cov), curvature = -0.5, linewidth = 0.2, lineend = "round") +
                scale_color_gradient(trans = "log10", low = "blue", high = "brown1", labels = label_number(accuracy = 1)) +
                scale_alpha(trans = "log10", range = c(0.01, 1)) + guides(alpha = "none") +
                coord_cartesian(xlim = c(0, max(exons$exon_end)), ylim = c(0, 3)) +
                theme_void() + theme(legend.position = "left")

        p2 <- ggplot(list_junctions[[i]], aes(x = donor, y = acceptor, color = avg_cov)) +
                theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                geom_vline(data = exons, aes(xintercept = exon_start), color = "lightgrey", linetype = "dashed") +
                geom_vline(data = exons, aes(xintercept = exon_end), color = "lightgrey", linetype = "dashed") +
                geom_hline(data = exons, aes(yintercept = exon_start), color = "lightgrey", linetype = "dashed") +
                geom_hline(data = exons, aes(yintercept = exon_end), color = "lightgrey", linetype = "dashed") +
                labs(x = "Start", y = "End", size = "Average") +
                geom_point(shape = 19, alpha = 1, size = 0.5) +
                scale_color_continuous(trans = "log10", low = "blue", high = "brown1", labels = label_number(accuracy = 1)) +
                scale_x_continuous(limits = c(-10, max(exons$exon_end)), breaks = seq(0, max(exons$exon_end), by = 100)) +  
                scale_y_continuous(limits = c(-10, max(exons$exon_end)), breaks = seq(0, max(exons$exon_end), by = 100)) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "left") +
                geom_rect(data = exons, inherit.aes = FALSE, fill = "darkgreen", alpha = 0.6, aes(xmin = exon_start, xmax = exon_end, ymin = -10, ymax = 0)) + 
                geom_rect(data = exons, inherit.aes = FALSE, fill = "darkgreen", alpha = 0.6, aes(ymin = exon_start, ymax = exon_end, xmin = -10, xmax = 0)) +
                geom_segment(data = introns, inherit.aes = FALSE, fill = "grey", alpha = 0.6, aes(x = intron_start, xend = intron_end, y = -5, yend = -5)) +
                geom_segment(data = introns, inherit.aes = FALSE, fill = "grey", alpha = 0.6, aes(y = intron_start, yend = intron_end, x = -5, xend = -5))

        p3 <- ggMarginal(p2, type = "density", margins = "both", color = "brown1", fill = "pink", size = 8)

        list_diagramplots[[i]] <- p1
        list_scatterplots[[i]] <- p3
    }

    return(list(list_junctions, list_diagramplots, list_scatterplots))
}
