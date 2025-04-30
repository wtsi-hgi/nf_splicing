workflow summarise_results_se {
    take:
    ch_sample

    main:


    emit:

}

process classify_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bed), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.txt"), path("${sample_id}.classified_junctions.reduce.txt") emit: ch_junctions
    tuple val(sample_id), path("${sample_id}.classified_junctions.png"), path("${sample_id}.classified_variants.png") emit: ch_pngs

    script:
    """
    ${projectDir}/scripts/classify_junctions.R -b ${bed} -e ${exon_pos} -p ${sample_id}
    """
}

