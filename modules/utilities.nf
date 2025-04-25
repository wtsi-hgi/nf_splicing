process cat_bed_file {
    label 'process_single'

    input:
    tuple val(sample_id), path(se_bed), path(pe_bed)

    output:
    tuple val(sample_id), path("${sample_id}.junctions.bed"), emit: ch_bed

    script:
    """
    cat ${se_bed} ${pe_bed} > ${sample_id}.junctions.bed
    """
}