workflow generate_summary_report {
    take:
    ch_sample

    main:
    ch_input_psi = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, idxstats, summary, canonical_barcodes, novel_barcodes, junctions, splicing_counts -> 
                                    tuple(sample, sample_id, splicing_counts) }
                            .groupTuple()

    ch_input_psi = ch_input_psi.filter { sample, sample_id, splicing_counts -> sample_id.size() == 3 }
                               .ifEmpty { exit 1, "Skipping creating html report due to the number of replicates after filtering." }

    CALCULATE_PSI(ch_input_psi)
    ch_psi_can_results = CALCULATE_PSI.out.ch_psi_can_results
    ch_psi_all_results = CALCULATE_PSI.out.ch_psi_all_results

    ch_report = ch_sample.map { sample_id, sample, barcode, exon_pos, merge_stats, trim_stats, idxstats, summary, canonical_barcodes, novel_barcodes, junctions, splicing_counts -> 
                                    tuple(sample, sample_id, barcode, exon_pos, merge_stats, trim_stats, idxstats, summary, canonical_barcodes, novel_barcodes, junctions) }
                         .groupTuple()
                         .join(ch_psi_can_results)
                         .join(ch_psi_all_results)

    CREATE_HTML_REPORT(ch_report)
    ch_junctions_category = CREATE_HTML_REPORT.out.ch_junctions_category
    ch_html_report = CREATE_HTML_REPORT.out.ch_html_report

    emit:
    ch_psi_can_results
    ch_psi_all_results
    ch_junctions_category
    ch_html_report
}

process CALCULATE_PSI {
    label 'process_single_dynamic_memory'
    
    memory {
        def file_size = splicing_counts[0].size()
        def mem = file_size <= 100_000_000 ? 4 :
                  file_size <= 1_000_000_000 ? 8 :
                  file_size <= 2_000_000_000 ? 16 :
                  file_size <= 4_000_000_000 ? 32 : 64
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/splicing_reports/${sample}", mode: "copy", overwrite: true

    input:
    tuple val(sample), val(sample_id), val(splicing_counts)

    output:
    tuple val(sample), path("${sample}.canon_only.corrected_psi.tsv"), emit: ch_psi_can_results
    tuple val(sample), path("${sample}.all_events.corrected_psi.tsv"), emit: ch_psi_all_results

    script:
    def list_sample_ids = sample_id.join(',')
    def list_splicing_counts = splicing_counts.join(',')

    """
    ${projectDir}/scripts/calculate_psi_with_error_model.R -r ${projectDir}/scripts \
                                                           -s ${list_sample_ids} \
                                                           -d ${list_splicing_counts} \
                                                           -c can \
                                                           -p ${sample}

    ${projectDir}/scripts/calculate_psi_with_error_model.R -r ${projectDir}/scripts \
                                                           -s ${list_sample_ids} \
                                                           -d ${list_splicing_counts} \
                                                           -c all \
                                                           -p ${sample}
    """
}

process CREATE_HTML_REPORT {
    label 'process_single_dynamic_memory'
    
    memory {
        def file_size_1 = canonical_barcodes[0].size()
        def file_size_2 = novel_barcodes[0].size()
        def file_size_total = file_size_1 + file_size_2
        def mem = file_size_total <= 100_000_000 ? 8 :
                  file_size_total <= 1_000_000_000 ? 16 :
                  file_size_total <= 2_000_000_000 ? 32 :
                  file_size_total <= 4_000_000_000 ? 64 : 128
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/splicing_reports/${sample}", mode: "copy", overwrite: true

    input:
    tuple val(sample), val(sample_id), val(barcode), val(exon_pos), 
          val(merge_stats), val(trim_stats), val(idxstats), val(summary), 
          val(canonical_barcodes), val(novel_barcodes), val(junctions), 
          val(psi_can_results), val(psi_all_results)

    output:
    tuple val(sample), path("${sample}.junctions_category.tsv"), emit: ch_junctions_category
    tuple val(sample), path("${sample}.splicing_report.html"), emit: ch_html_report

    script:
    def list_sample_ids = sample_id.join(',')
    def file_barcode = barcode[0]
    def file_exon_pos = exon_pos[0]
    def list_merge_stats = merge_stats.join(',')
    def list_trim_stats = trim_stats.join(',')
    def list_idxstats = idxstats.join(',')
    def list_summary = summary.join(',')
    def list_canonical_barcodes = canonical_barcodes.join(',')
    def list_novel_barcodes = novel_barcodes.join(',')
    def list_junctions = junctions.join(',')
    def file_psi_can_results = psi_can_results
    def file_psi_all_results = psi_all_results

    """
    ${projectDir}/scripts/create_html_report.R -r ${projectDir}/scripts \
                                               -l ${params.library} \
                                               -b ${file_barcode} \
                                               -e ${file_exon_pos} \
                                               -s ${list_sample_ids} \
                                               -t ${list_trim_stats} \
                                               -m ${list_merge_stats} \
                                               -f ${list_idxstats} \
                                               -a ${list_summary} \
                                               -c ${list_canonical_barcodes} \
                                               -n ${list_novel_barcodes} \
                                               -j ${list_junctions} \
                                               -d ${file_psi_can_results},${file_psi_all_results} \
                                               -p ${sample}
    """
}