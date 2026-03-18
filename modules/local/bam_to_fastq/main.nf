process SE_BAM_TO_FASTQ {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.wrongmap.fastq.gz"), emit: ch_se_wrongmap

    script:
    """
    samtools fastq -@ 32 -c 9 -0 ${sample_id}.bwa_se.wrongmap.fastq.gz ${bam}
    rm ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}

process PE_BAM_TO_FASTQ {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.wrongmap.r1.fastq.gz"), path("${sample_id}.bwa_pe.wrongmap.r2.fastq.gz"), emit: ch_pe_wrongmap

    script:
    """
    samtools fastq -@ 32 -c 9 -1 ${sample_id}.bwa_pe.wrongmap.r1.fastq.gz -2 ${sample_id}.bwa_pe.wrongmap.r2.fastq.gz -n ${bam}
    rm ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}
