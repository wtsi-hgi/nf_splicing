workflow summarise_results {
    take:
    ch_sample

    main:


    emit:

}

process cat_bed_file {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(se_bed), path(pe_bed)

    output:
    tuple val(sample_id), path("${sample_id}.junctions.bed"), emit: ch_bed

    script:
    """
    cat ${se_bed} ${pe_bed} > ${sample_id}.junctions.bed
    """
}

process cat_spliced_file {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(se_spliced), path(pe_spliced)

    output:
    tuple val(sample_id), path("${sample_id}.spliced_products.txt"), emit: ch_spliced

    script:
    """
    cat ${se_spliced} ${pe_spliced} > ${sample_id}.spliced_products.txt
    """
}

process classify_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bed), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.junctions.txt"), path("${sample_id}.junctions.reduce.txt"), emit: ch_classified_junctions
    tuple val(sample_id), path("${sample_id}.junctions.png"), path("${sample_id}.variants.png"), emit: ch_classified_png

    script:
    """
    ${projectDir}/scripts/classify_junctions.R -b ${bed} -e ${exon_pos} -m ${params.classify_min_overlap} -c ${params.classify_min_cov} -r ${params.classify_reduce}
    """
}