create_barplots <- function(sample_reps, total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads)
{
    summary_reads <- as.data.table(cbind(sample_reps, total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads))
    colnames(summary_reads) <- c("sample_reps", "total_reads", "merged_reads", "unmerged_reads", "inclusion_reads", "skipping_reads", "map_reads", "unexplain_reads")

    summary_reads <- summary_reads %>% mutate(across(c(total_reads, merged_reads, unmerged_reads, inclusion_reads, skipping_reads, map_reads, unexplain_reads), as.numeric))
    summary_pct <- summary_reads %>%
                    mutate(pct_merged = round((merged_reads / total_reads) * 100, 1),
                           pct_unmerged = round((unmerged_reads / total_reads) * 100, 1),
                           library_coverage = as.integer(total_reads / length(unique(barcode_association$varid))),
                           pct_inclusion = round((inclusion_reads / total_reads) * 100, 1),
                           pct_skipping = round((skipping_reads / total_reads) * 100, 1),
                           pct_map = round((map_reads / total_reads) * 100, 1),
                           pct_unexplain = round((unexplain_reads / total_reads) * 100, 1))

    tib_pct_total <- summary_pct %>% 
                        select(sample_reps, pct_merged, pct_unmerged) %>%
                        pivot_longer(cols = -sample_reps, names_to = "type", values_to = "pct")
    tib_pct_total$type <- factor(tib_pct_total$type, levels = c("pct_unmerged", "pct_merged"))
    tib_cov <- summary_reads %>% select(sample_reps, library_coverage)
    tib_cov$type <- "coverage"

    tib_pct_canonical <- summary_pct %>% 
                            select(sample_reps, pct_inclusion, pct_skipping) %>%
                            pivot_longer(cols = -sample_reps, names_to = "type", values_to = "pct")

    tib_pct_novel <- summary_reads %>% 
                        select(sample_reps, pct_map, pct_unexplain) %>%
                        pivot_longer(cols = -sample_reps, names_to = "type", values_to = "pct")

    select_colors <- select_colorblind("col15")[1:2]
    fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

    p1 <- ggplot(tib_pct_total,  aes(x = reps, y = pct, fill = type)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_fill_manual(values = fill_colors) +
            labs(x = NULL, y = "percent") +
            theme(legend.position = "right", legend.title = element_blank()) +
            theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
            theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
            theme(axis.text = element_text(size = 8, face = "bold")) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)

    p2 <- ggplot(tib_pct_canonical,  aes(x = reps, y = pct, fill = type)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_fill_manual(values = fill_colors) +
            labs(x = NULL, y = "percent") +
            theme(legend.position = "right", legend.title = element_blank()) +
            theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
            theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
            theme(axis.text = element_text(size = 8, face = "bold")) +
            theme(axis.text.x = element_text(angle = 90)) +
            geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3)

    p3 <- ggplot(tib_pct_novel,  aes(x = reps, y = pct, fill = type)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_fill_manual(values = fill_colors) +
            labs(x = NULL, y = "percent", title = sample_prefix) +
            theme(legend.position = "right", legend.title = element_blank()) +
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
