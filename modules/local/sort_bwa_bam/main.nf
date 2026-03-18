process SORT_BWA_SE_BAM {
    label 'process_high_memory'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.filtered.sorted.bam"), path("${sample_id}.bwa_se.filtered.sorted.bam.bai"), emit: ch_se_sorted_bam
    tuple val(sample_id), path("${sample_id}.bwa_se.canonical_stats.tsv"), emit: ch_se_canonical_stats

    script:
    """
    samtools sort -@ 40 -o ${sample_id}.bwa_se.filtered.sorted.bam ${bam}
    samtools index -@ 40 ${sample_id}.bwa_se.filtered.sorted.bam
    samtools idxstats ${sample_id}.bwa_se.filtered.sorted.bam > ${sample_id}.bwa_se.canonical_stats.tsv
    rm ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}

process SORT_BWA_PE_BAM {
    label 'process_high_memory'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.filtered.sorted.bam"), path("${sample_id}.bwa_pe.filtered.sorted.bam.bai"), emit: ch_pe_sorted_bam
    tuple val(sample_id), path("${sample_id}.bwa_pe.canonical_stats.tsv"), emit: ch_pe_canonical_stats

    script:
    """
    samtools sort -@ 40 -o ${sample_id}.bwa_pe.filtered.sorted.bam ${bam}
    samtools index -@ 40 ${sample_id}.bwa_pe.filtered.sorted.bam
    samtools idxstats ${sample_id}.bwa_pe.filtered.sorted.bam > ${sample_id}.bwa_pe.canonical_stats.tsv
    rm ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}
