process SORT_SE_BAM {
    label 'process_high_memory'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.fixed.sorted.bam"), path("${sample_id}.hisat2_se.fixed.sorted.bam.bai"), emit: ch_se_sorted_bam 

    script:
    """
    samtools sort -@ 40 -o ${sample_id}.hisat2_se.fixed.sorted.bam ${bam}
    samtools index -@ 40 ${sample_id}.hisat2_se.fixed.sorted.bam
    bamtools stats -in ${sample_id}.hisat2_se.fixed.sorted.bam > ${sample_id}.hisat2_se.fixed.tsv
    rm ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}

process SORT_PE_BAM {
    label 'process_high_memory'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.fixed.sorted.bam"), path("${sample_id}.hisat2_pe.fixed.sorted.bam.bai"), emit: ch_pe_sorted_bam 

    script:
    """
    samtools sort -@ 40 -o ${sample_id}.hisat2_pe.fixed.sorted.bam ${bam}
    samtools index -@ 40 ${sample_id}.hisat2_pe.fixed.sorted.bam
    bamtools stats -in ${sample_id}.hisat2_pe.fixed.sorted.bam > ${sample_id}.hisat2_pe.fixed.tsv
    rm ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}
