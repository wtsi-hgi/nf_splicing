workflow summarise_results {
    take:
    ch_sample

    main:
    /* -- 1. classify junctions -- */
    ch_junctions = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                    tuple(sample_id, exon_pos, junctions) }
    classify_junctions(ch_junctions)
    ch_classified_junctions = classify_junctions.out.ch_classified_junctions
    ch_classified_plots = classify_junctions.out.ch_classified_plots

    /* -- 2. create junction plots -- */
    ch_plots = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                tuple(sample_id, sample) }
                        .join(ch_classified_junctions.map { sample_id, classified_junctions, reduced_junctions -> tuple(sample_id, classified_junctions) })
                        .map { sample_id, sample, classified_junctions ->  tuple(sample, sample_id, classified_junctions)}
                        .groupTuple()
    ch_exon_pos = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                tuple(sample, exon_pos) }
                           .distinct()
    ch_plots = ch_plots.join(ch_exon_pos)

    create_junction_plots(ch_plots)
    ch_junction_plots = create_junction_plots.out.ch_junction_plots

    /* -- 3. create count matrix -- */
    ch_input = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                tuple(sample_id, barcode, filter_barcodes) }
                        .join(ch_classified_junctions.map { sample_id, classified_junctions, reduced_junctions -> tuple(sample_id, classified_junctions) })
    create_splicing_matrix(ch_input)
    ch_splicing_matrix = create_splicing_matrix.out.ch_splicing_matrix

    /* -- 4. create summary report -- */
    ch_barcode_association = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                                tuple(sample, barcode) }
                                      .distinct()
    ch_report = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                    tuple(sample, sample_id, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary) }
                         .groupTuple()
                         
    // Check if the number of replicates equals 3, otherwise skip the html report
    ch_report_filtered = ch_report.filter { sample, sample_id, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary ->
                                                merge_stats.size() == 3 }
                                  .ifEmpty { exit 1, "Skipping creating html report due to the number of replicates." }

    ch_report_filtered = ch_barcode_association.join(ch_report_filtered)
                                               .join(ch_junction_plots)
                                               .map{ sample, barcode, sample_id, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, 
                                                        junction_venn, junction_scatter, junction_view, junction_distribution, junction_corr ->
                                                     tuple(sample, barcode, sample_id, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, 
                                                        [junction_venn, junction_scatter, junction_view, junction_distribution, junction_corr]) }
    
    create_html_report(ch_report_filtered)
    ch_html_report = create_html_report.out.ch_html_report

    emit:
    ch_classified_junctions
    ch_classified_plots
    ch_junction_plots
    ch_splicing_matrix
    ch_html_report
}

process classify_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_junctions/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(exon_pos), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.txt"), path("${sample_id}.classified_junctions.reduce.txt"), emit: ch_classified_junctions
    tuple val(sample_id), path("${sample_id}.classified_junctions.png"), path("${sample_id}.classified_variants.png"), emit: ch_classified_plots

    script:
    """
    ${projectDir}/scripts/classify_junctions.R -b ${bed} -e ${exon_pos} -m ${params.classify_min_overlap} -c ${params.classify_min_cov} -r ${params.classify_reduce} -p ${sample_id}
    """
}

process create_junction_plots {
    label 'process_single'

    publishDir "${params.outdir}/splicing_reports/${sample}", mode: "copy", overwrite: true

    input:
    tuple val(sample), val(sample_id), val(classified_junctions), path(exon_pos)

    output:
    tuple val(sample), path("${sample}.junction_venn.png"), 
                       path("${sample}.junction_scatter.png"), 
                       path("${sample}.junction_view.png"),
                       path("${sample}.junction_distribution.png"), 
                       path("${sample}.junction_corr.png"), emit: ch_junction_plots

    script:
    """
    ${projectDir}/scripts/create_junction_plots.R -s ${sample_id.join(',')} -n ${classified_junctions.join(',')} -e ${exon_pos} -p ${sample}
    """
}

process create_splicing_matrix {
    label 'process_single'

    publishDir "${params.outdir}/splicing_counts", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), path(filter_barcodes), path(classified_junctions)

    output:
    tuple val(sample_id), path("${sample_id}.splicing_matrix.txt"), emit: ch_splicing_matrix

    script:
    """
    ${projectDir}/scripts/create_splicing_matrix.R -b ${barcode} -c ${filter_barcodes} -n ${classified_junctions} -p ${sample_id}
    """
}

process create_html_report {
    label 'process_single'

    publishDir "${params.outdir}/splicing_reports/${sample}", mode: "copy", overwrite: true

    input:
    tuple val(sample), path(barcode), val(sample_id), 
          val(merge_stats), val(trim_stats), val(sample_idxstats), val(filter_barcodes), val(map_barcodes), val(summary), 
          val(junction_plots)

    output:
    tuple val(sample), path("${sample}.splicing_report.html"), emit: ch_html_report
    
    script:
    def list_sample_ids = sample_id.join(',')
    def list_trim_stats = trim_stats.join(',')
    def list_merge_stats = merge_stats.join(',')
    def list_sample_idxstats = sample_idxstats.join(',')
    def list_filter_barcodes = filter_barcodes.join(',')
    def list_map_barcodes = map_barcodes.join(',')
    def list_summary = summary.join(',')
    def list_junction_plots = junction_plots.join(',')

    """
    ${projectDir}/scripts/create_html_report.R -b ${barcode} \
                                               -s ${list_sample_ids} \
                                               -t ${list_trim_stats} \
                                               -m ${list_merge_stats} \
                                               -f ${list_sample_idxstats} \
                                               -a ${list_summary} \
                                               -c ${list_filter_barcodes} \
                                               -n ${list_map_barcodes} \
                                               -j ${list_junction_plots} \
                                               -p ${sample}
    """
}