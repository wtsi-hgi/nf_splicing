process BWA_ALIGN_SE_READS {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(extended_frags), path(exon_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.unique.bam"), emit: ch_se_unique_bam
    tuple val(sample_id), path("${sample_id}.bwa_se.unmapped.fastq.gz"), emit: ch_se_unmapped

    script:
    """
    bwa index ${exon_fasta}
    bwa mem -t 32 -B ${params.bwa_mismatch} \
                  -O ${params.bwa_gap_open} \
                  -E ${params.bwa_gap_ext} \
                  -L ${params.bwa_clip} \
                  ${exon_fasta} ${extended_frags} > ${sample_id}.bwa_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.bwa_se.unmapped.fastq.gz ${sample_id}.bwa_se.sam
    samtools view -@ 32 -b -F 4 -F 256 -F 2048 ${sample_id}.bwa_se.sam > ${sample_id}.bwa_se.unique.bam
    rm ${sample_id}.bwa_se.sam
    """
}

process BWA_ALIGN_PE_READS {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(not_combined_1), path(not_combined_2), path(exon_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.unique.bam"), emit: ch_pe_unique_bam
    tuple val(sample_id), path("${sample_id}.bwa_pe.unmapped.r1.fastq.gz"), path("${sample_id}.bwa_pe.unmapped.r2.fastq.gz"), emit: ch_pe_unmapped

    script:
    """
    bwa index ${exon_fasta}
    bwa mem -t 32 -O ${params.bwa_gap_open} \
                  -E ${params.bwa_gap_ext} \
                  -L ${params.bwa_clip} \
                  ${exon_fasta} ${not_combined_1} ${not_combined_2} > ${sample_id}.bwa_pe.sam
    samtools fastq -@ 32 -F 2 -c 9 -1 ${sample_id}.bwa_pe.unmapped.r1.fastq.gz -2 ${sample_id}.bwa_pe.unmapped.r2.fastq.gz -n ${sample_id}.bwa_pe.sam
    samtools view -b -f 2 -F 256 -F 2048 ${sample_id}.bwa_pe.sam > ${sample_id}.bwa_pe.unique.bam
    rm ${sample_id}.bwa_pe.sam
    """
}
