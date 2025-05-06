workflow summarise_results {
    take:
    ch_sample

    main:
    ch_junctions = ch_sample.map { sample_id, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                    tuple(sample_id, exon_pos, junctions) }
    classify_junctions(ch_junctions)
    ch_classified_junctions = classify_junctions.out.ch_classified_junctions
    ch_classified_png = classify_junctions.out.ch_classified_png

    ch_input = ch_sample.map { sample_id, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                tuple(sample_id, barcode, filter_barcodes) }
                        .join(ch_classified_junctions.map { sample_id, classified_junctions, reduced_junctions -> tuple(sample_id, classified_junctions) })
    create_splicing_matrix(ch_input)
    ch_splicing_matrix = create_splicing_matrix.out.ch_splicing_matrix

    ch_report = ch_sample.map { sample_id, barcode, exon_pos, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary, junctions -> 
                                    tuple(sample_id, barcode, merge_stats, trim_stats, sample_idxstats, filter_barcodes, map_barcodes, summary) }
    ch_report.collect().view()


    emit:
    ch_classified_junctions
    ch_classified_png
    ch_splicing_matrix
}

process classify_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(exon_pos), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.txt"), path("${sample_id}.classified_junctions.reduce.txt"), emit: ch_classified_junctions
    tuple val(sample_id), path("${sample_id}.classified_junctions.png"), path("${sample_id}.classified_variants.png"), emit: ch_classified_png

    script:
    """
    ${projectDir}/scripts/classify_junctions.R -b ${bed} -e ${exon_pos} -m ${params.classify_min_overlap} -c ${params.classify_min_cov} -r ${params.classify_reduce} -p ${sample_id}
    """
}

process create_splicing_matrix {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), path(filter_barcodes), path(classified_junctions)

    output:
    tuple val(sample_id), path("${sample_id}.splicing_matrix.txt"), emit: ch_splicing_matrix

    script:
    """
    ${projectDir}/scripts/create_splicing_matrix.R -b ${barcode} -c ${filter_barcodes} -n ${classified_junctions} -p ${sample_id}
    """
}