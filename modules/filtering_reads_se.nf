workflow filtering_reads {
    take:
    ch_sample

    main:


    emit:
}

process bwa_align_se_reads {
    label 'process_low'

    input:
    tuple val(sample_id), path(barcode), path(extended_frags), path(exon_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.fail.fastq.gz"), emit: ch_bwa_se

    script:
    """
    bwa mem -t 32 -O 10,10 -E 5,5 -L 1,1 ${exon_fasta} ${extended_frags} > ${sample_id}.filter_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.filter_se.unmapped.fastq.gz ${sample_id}.filter_se.sam
    samtools view -@ 32 -b -F 4 -F 256 -F 2048 ${sample_id}.filter_se.sam > ${sample_id}.filter_se.unique.bam
    samtools sort -@ 32 -o ${sample_id}.filter_se.unique.sorted.bam ${sample_id}.filter_se.unique.bam
    samtools index ${sample_id}.filter_se.unique.sorted.bam
    rm ${sample_id}.filter_se.sam ${sample_id}.filter_se.unique.bam

    """
}

