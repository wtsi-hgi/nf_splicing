workflow filtering_reads_se {
    take:
    ch_sample

    main:
    bwa_align_se_reads(ch_sample)
    ch_bwa_se_unmapped = bwa_align_se_reads.out.ch_bwa_se_unmapped
    ch_bwa_se_bam = bwa_align_se_reads.out.ch_bwa_se_bam
    ch_exon_pos = bwa_align_se_reads.out.ch_exon_pos
    
    ch_bwa_se_bam = ch_bwa_se_bam.join(ch_exon_pos)
    filter_se_reads(ch_bwa_se_bam)
    ch_bwa_se_wrongmap = filter_se_reads.out.ch_bwa_se_wrongmap
    ch_bwa_se_filtered = filter_se_reads.out.ch_bwa_se_filtered
    ch_bwa_se_filtered_idxstats = filter_se_reads.out.ch_bwa_se_filtered_idxstats

    ch_bwa_se_joined = ch_bwa_se_unmapped.join(ch_bwa_se_wrongmap)
    merge_fastqs(ch_bwa_se_joined)
    ch_bwa_se_fail = merge_fastqs.out.ch_bwa_se_fail

    extract_barcodes_se(ch_bwa_se_filtered)
    ch_bwa_se_barcodes = extract_barcodes_se.out.ch_bwa_se_barcodes

    emit:
    ch_exon_pos
    ch_bwa_se_filtered
    ch_bwa_se_filtered_idxstats
    ch_bwa_se_fail
    ch_bwa_se_barcodes
}

process bwa_align_se_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(barcode), path(extended_frags), path(exon_fasta), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.unmapped.fastq.gz"), emit: ch_bwa_se_unmapped
    tuple val(sample_id), path("${sample_id}.filter_se.unique.sorted.bam"), path("${sample_id}.filter_se.unique.sorted.bam.bai"), emit: ch_bwa_se_bam
    tuple val(sample_id), path(exon_pos), emit: ch_exon_pos

    script:
    """
    bwa index ${exon_fasta}
    bwa mem -t 32 -O ${params.bwa_gap_open} -E ${params.bwa_gap_ext} -L ${params.bwa_clip} ${exon_fasta} ${extended_frags} > ${sample_id}.filter_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.filter_se.unmapped.fastq.gz ${sample_id}.filter_se.sam
    samtools view -@ 32 -b -F 4 -F 256 -F 2048 ${sample_id}.filter_se.sam > ${sample_id}.filter_se.unique.bam
    samtools sort -@ 32 -o ${sample_id}.filter_se.unique.sorted.bam ${sample_id}.filter_se.unique.bam
    samtools index -@ 32 ${sample_id}.filter_se.unique.sorted.bam
    rm ${sample_id}.filter_se.sam ${sample_id}.filter_se.unique.bam
    """
}

process filter_se_reads {
    label 'process_high'

    input:
    tuple val(sample_id), path(bam), path(bai), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.wrongmap.fastq.gz"), emit: ch_bwa_se_wrongmap
    tuple val(sample_id), path("${sample_id}.filter_se.filtered.sorted.bam"), path("${sample_id}.filter_se.filtered.sorted.bam.bai"), emit: ch_bwa_se_filtered
    tuple val(sample_id), path("${sample_id}.filter_se.filtered.idxstats.txt"), emit: ch_bwa_se_filtered_idxstats

    script:
    """
    ${projectDir}/scripts/filter_bam_se.R -b ${bam} -e ${exon_pos}
    samtools view -@ 80 -b -o ${sample_id}.filter_se.filtered.bam ${sample_id}.filter_se.unique.sorted.filtered.sam
    samtools sort -@ 80 -o ${sample_id}.filter_se.filtered.sorted.bam ${sample_id}.filter_se.filtered.bam
    samtools index -@ 80 ${sample_id}.filter_se.filtered.sorted.bam
    rm ${sample_id}.filter_se.filtered.bam ${sample_id}.filter_se.unique.sorted.filtered.sam
    samtools idxstats ${sample_id}.filter_se.filtered.sorted.bam > ${sample_id}.filter_se.filtered.idxstats.txt

    mv ${sample_id}.filter_se.unique.sorted.wrongmap.fastq ${sample_id}.filter_se.wrongmap.fastq
    pigz --best -p 80 ${sample_id}.filter_se.wrongmap.fastq
    """
}

process merge_fastqs {
    label 'process_single'

    input:
    tuple val(sample_id), path(unmapped_fastq), path(wrongmap_fastq)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.fail.fastq.gz"), emit: ch_bwa_se_fail

    script:
    """
    cat ${unmapped_fastq} ${wrongmap_fastq} > ${sample_id}.filter_se.fail.fastq.gz
    """
}

process extract_barcodes_se {
    label 'process_medium'

    publishDir "${params.outdir}/extracted_barcodes", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.barcodes.txt"), emit: ch_bwa_se_barcodes

    script:
    """
    ${projectDir}/scripts/extract_barcodes_se.R -b ${bam} -p ${params.barcode_template} -m ${params.barcode_marker}
    mv ${sample_id}.filter_se.filtered.sorted.barcodes.txt ${sample_id}.filter_se.barcodes.txt
    """
}