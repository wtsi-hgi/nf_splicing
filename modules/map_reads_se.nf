workflow map_reads_se {
    take:
    ch_sample

    main:
    hisat2_align_se_reads(ch_sample)
    ch_hisat2_se_unmapped = hisat2_align_se_reads.out.ch_hisat2_se_unmapped
    ch_hisat2_se_bam = hisat2_align_se_reads.out.ch_hisat2_se_bam
    ch_exon_pos = hisat2_align_se_reads.out.ch_exon_pos

    ch_hisat2_se_bam_join = ch_sample.map { sample_id, barcode, exon_pos, read, ref_fasta ->
                                            tuple(sample_id, barcode, ref_fasta) }
                                     .join(ch_hisat2_se_bam)

    emit:
}

process hisat2_align_se_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(barcode), path(exon_pos), path(read), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.unmapped.fastq.gz"), emit: ch_hisat2_se_unmapped
    tuple val(sample_id), path("${sample_id}.map_se.unique.sorted.bam"), path("${sample_id}.map_se.unique.sorted.bam.bai"), emit: ch_hisat2_se_bam
    tuple val(sample_id), path(exon_pos), emit: ch_exon_pos

    script:
    """
    hisat2-build ${ref_fasta} ${ref_fasta}
    hisat2 -x ${ref_fasta} \
           -U ${read} \
           --score-min ${params.hisat2_score_min} \
           --mp ${params.hisat2_mp} \
           --sp ${params.hisat2_sp} \
           --np ${params.hisat2_np} \
           --pen-noncansplice ${params.hisat2_pen_noncansplice} \
           --summary-file ${sample_id}.map_se.summary.txt \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.map_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.map_se.unmapped.fastq.gz ${sample_id}.map_se.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==0)||(\$2==16)){print \$0}}' ${sample_id}.map_se.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.map_se.unique.bam
    samtools sort -@ 32 -o ${sample_id}.map_se.unique.sorted.bam ${sample_id}.map_se.unique.bam
    samtools index -@ 32 ${sample_id}.map_se.unique.sorted.bam
    bamtools stats -in ${sample_id}.map_se.unique.sorted.bam > ${sample_id}.map_se.unique.txt
    rm ${sample_id}.map_se.sam ${sample_id}.map_se.unique.bam
    """
}

process fix_se_reads {
    label 'process_high'

    input:
    

    output:


    script:
    """

    """
}
