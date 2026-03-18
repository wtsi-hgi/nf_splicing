process HISAT2_ALIGN_SE_READS {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(se_reads), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.unique.bam"), emit: ch_se_unique_bam
    tuple val(sample_id), path("${sample_id}.hisat2_se.unmapped.fastq.gz"), emit: ch_se_unmapped
    tuple val(sample_id), path("${sample_id}.hisat2_se.novel_stats.tsv"), emit: ch_se_novel_stats

    script:
    """
    hisat2-build ${ref_fasta} ${ref_fasta}
    hisat2 -x ${ref_fasta} \
           -U ${se_reads} \
           --score-min ${params.hisat2_score_min} \
           --mp ${params.hisat2_mp} \
           --sp ${params.hisat2_sp} \
           --np ${params.hisat2_np} \
           --pen-noncansplice ${params.hisat2_pen_noncansplice} \
           --summary-file ${sample_id}.hisat2_se.novel_stats.tsv \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.hisat2_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.hisat2_se.unmapped.fastq.gz ${sample_id}.hisat2_se.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==0)||(\$2==16)){print \$0}}' ${sample_id}.hisat2_se.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.hisat2_se.unique.bam
    rm ${sample_id}.hisat2_se.sam
    """
}

process HISAT2_ALIGN_PE_READS {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.unique.bam"), emit: ch_pe_unique_bam
    tuple val(sample_id), path("${sample_id}.hisat2_pe.unmapped.R1.fastq.gz"), path("${sample_id}.hisat2_pe.unmapped.R2.fastq.gz"), emit: ch_pe_unmapped
    tuple val(sample_id), path("${sample_id}.hisat2_pe.novel_stats.tsv"), emit: ch_pe_novel_stats

    script:
    """
    hisat2-build ${ref_fasta} ${ref_fasta}
    hisat2 -x ${ref_fasta} \
           -1 ${read1} -2 ${read2} --fr \
           --score-min ${params.hisat2_score_min} \
           --mp ${params.hisat2_mp} \
           --sp ${params.hisat2_sp} \
           --np ${params.hisat2_np} \
           --pen-noncansplice ${params.hisat2_pen_noncansplice} \
           --summary-file ${sample_id}.hisat2_pe.novel_stats.tsv \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.hisat2_pe.sam
    samtools fastq -@ 32 -F 2 -c 9 -1 ${sample_id}.hisat2_pe.unmapped.R1.fastq.gz -2 ${sample_id}.hisat2_pe.unmapped.R2.fastq.gz -n ${sample_id}.hisat2_pe.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==99)||(\$2==147)||(\$2==83)||(\$2==163)){print \$0}}' ${sample_id}.hisat2_pe.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.hisat2_pe.unique.bam
    rm ${sample_id}.hisat2_pe.sam
    """
}
