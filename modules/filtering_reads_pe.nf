workflow filtering_reads_pe {
    take:
    ch_sample

    main:
    bwa_align_pe_reads(ch_sample)
    ch_bwa_pe_unmapped = bwa_align_pe_reads.out.ch_bwa_pe_unmapped
    ch_bwa_pe_bam = bwa_align_pe_reads.out.ch_bwa_pe_bam
    ch_exon_pos = bwa_align_pe_reads.out.ch_exon_pos
    
    ch_bwa_pe_bam = ch_bwa_pe_bam.join(ch_exon_pos)
    filter_pe_reads(ch_bwa_pe_bam)
    ch_bwa_pe_wrongmap = filter_pe_reads.out.ch_bwa_pe_wrongmap
    ch_bwa_pe_filtered = filter_pe_reads.out.ch_bwa_pe_filtered

    ch_bwa_pe_joined = ch_bwa_pe_unmapped.join(ch_bwa_pe_wrongmap)
    merge_fastqs(ch_bwa_pe_joined)
    ch_bwa_pe_fail = merge_fastqs.out.ch_bwa_pe_fail

    emit:
    ch_bwa_pe_filtered
    ch_bwa_pe_fail
    ch_exon_pos
}

process bwa_align_pe_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(barcode), path(not_combined_1), path(not_combined_2), path(exon_fasta), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.filter_pe.unmapped.r1.fastq.gz"), path("${sample_id}.filter_pe.unmapped.r2.fastq.gz"), emit: ch_bwa_pe_unmapped
    tuple val(sample_id), path("${sample_id}.filter_pe.unique.sorted.bam"), path("${sample_id}.filter_pe.unique.sorted.bam.bai"), emit: ch_bwa_pe_bam
    tuple val(sample_id), path(exon_pos), emit: ch_exon_pos

    script:
    """
    bwa index ${exon_fasta}
    bwa mem -t 32 -O 10,10 -E 5,5 -L 1,1 ${exon_fasta} ${not_combined_1} ${not_combined_2} > ${sample_id}.filter_pe.sam
    samtools fastq -@ 32 -F 2 -c 9 -1 ${sample_id}.filter_pe.unmapped.r1.fastq.gz -2 ${sample_id}.filter_pe.unmapped.r2.fastq.gz -n ${sample_id}.filter_pe.sam
    samtools view -b -f 2 -F 256 -F 2048 ${sample_id}.filter_pe.sam > ${sample_id}.filter_pe.unique.bam
    samtools sort -@ 32 -o ${sample_id}.filter_pe.unique.sorted.bam ${sample_id}.filter_pe.unique.bam
    samtools index -@ 32 ${sample_id}.filter_pe.unique.sorted.bam
    rm ${sample_id}.filter_pe.sam ${sample_id}.filter_pe.unique.bam
    """
}

process filter_pe_reads {
    label 'process_high'

    input:
    tuple val(sample_id), path(bam), path(bai), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.filter_pe.wrongmap.r1.fastq.gz"), path("${sample_id}.filter_pe.wrongmap.r2.fastq.gz"), emit: ch_bwa_pe_wrongmap
    tuple val(sample_id), path("${sample_id}.filter_pe.filtered.sorted.bam"), path("${sample_id}.filter_pe.filtered.sorted.bam.bai"), emit: ch_bwa_pe_filtered

    script:
    """
    ${projectDir}/scripts/filter_bam_pe.R -b ${bam} -e ${exon_pos}
    samtools view -@ 80 -b -o ${sample_id}.filter_pe.filtered.bam ${sample_id}.filter_pe.unique.sorted.filtered.sam
    samtools sort -@ 80 -o ${sample_id}.filter_pe.filtered.sorted.bam ${sample_id}.filter_pe.filtered.bam
    samtools index -@ 80 ${sample_id}.filter_pe.filtered.sorted.bam
    rm ${sample_id}.filter_pe.filtered.bam ${sample_id}.filter_pe.unique.sorted.filtered.sam
    
    mv ${sample_id}.filter_pe.unique.sorted.wrongmap.r1.fastq ${sample_id}.filter_pe.wrongmap.r1.fastq
    pigz --best -p 80 ${sample_id}.filter_pe.wrongmap.r1.fastq
    mv ${sample_id}.filter_pe.unique.sorted.wrongmap.r2.fastq ${sample_id}.filter_pe.wrongmap.r2.fastq
    pigz --best -p 80 ${sample_id}.filter_pe.wrongmap.r2.fastq
    """
}

process merge_fastqs {
    label 'process_single'

    input:
    tuple val(sample_id), path(unmapped_fastq_r1), path(unmapped_fastq_r2), path(wrongmap_fastq_r1), path(wrongmap_fastq_r2)

    output:
    tuple val(sample_id), path("${sample_id}.filter_pe.fail.r1.fastq.gz"), path("${sample_id}.filter_pe.fail.r2.fastq.gz"), emit: ch_bwa_pe_fail

    script:
    """
    cat ${unmapped_fastq_r1} ${wrongmap_fastq_r1} > ${sample_id}.filter_pe.fail.r1.fastq.gz
    cat ${unmapped_fastq_r2} ${wrongmap_fastq_r2} > ${sample_id}.filter_pe.fail.r2.fastq.gz
    """
}
