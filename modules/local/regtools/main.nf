process EXTRACT_SE_JUNCTIONS {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.junctions.bed"), emit: ch_se_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.hisat2_se.junctions.bed ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$( regtools --version 2>&1 | grep -i "Version" | awk '{print \$2}' )
    END_VERSIONS
    """
}

process EXTRACT_PE_JUNCTIONS {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.junctions.bed"), emit: ch_pe_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.hisat2_pe.junctions.bed ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$( regtools --version 2>&1 | grep -i "Version" | awk '{print \$2}' )
    END_VERSIONS
    """
}
