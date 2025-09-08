workflow filter_reads_pe {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    bwa_align_pe_reads(ch_sample)
    ch_bwa_pe_unmapped = bwa_align_pe_reads.out.ch_bwa_pe_unmapped
    ch_bwa_pe_bam = bwa_align_pe_reads.out.ch_bwa_pe_bam
    ch_exon_pos = bwa_align_pe_reads.out.ch_exon_pos
    
    /* -- 2. filter reads -- */
    ch_bwa_pe_bam = ch_bwa_pe_bam.join(ch_exon_pos)
    filter_pe_reads(ch_bwa_pe_bam)
    ch_bwa_pe_wrongmap = filter_pe_reads.out.ch_bwa_pe_wrongmap
    ch_bwa_pe_filtered = filter_pe_reads.out.ch_bwa_pe_filtered
    ch_bwa_pe_filtered_idxstats = filter_pe_reads.out.ch_bwa_pe_filtered_idxstats

    /* -- 3. cat reads -- */
    ch_bwa_pe_joined = ch_bwa_pe_unmapped.join(ch_bwa_pe_wrongmap)
    merge_fastqs(ch_bwa_pe_joined)
    ch_bwa_pe_fail = merge_fastqs.out.ch_bwa_pe_fail

    /* -- 4. extract barcodes -- */
    ch_bwa_pe_filtered_joined = ch_bwa_pe_filtered.join(ch_sample.map { sample_id, barcode, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                                        tuple(sample_id, barcode) })
    extract_barcodes_pe(ch_bwa_pe_filtered_joined)
    ch_bwa_pe_barcodes = extract_barcodes_pe.out.ch_bwa_pe_barcodes

    emit:
    ch_bwa_pe_filtered
    ch_bwa_pe_filtered_idxstats
    ch_bwa_pe_fail
    ch_bwa_pe_barcodes
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
    bwa mem -t 32 -O ${params.bwa_gap_open} \
                  -E ${params.bwa_gap_ext} \
                  -L ${params.bwa_clip} \
                  ${exon_fasta} ${not_combined_1} ${not_combined_2} > ${sample_id}.filter_pe.sam
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
    tuple val(sample_id), path("${sample_id}.filter_pe.filtered.idxstats.txt"), emit: ch_bwa_pe_filtered_idxstats

    script:
    """
    ${projectDir}/scripts/filter_bam_pe.R -b ${bam} -e ${exon_pos} -s ${params.filter_softclip_base}
    samtools view -@ 64 -b -o ${sample_id}.filter_pe.filtered.bam ${sample_id}.filter_pe.unique.sorted.filtered.sam
    samtools sort -@ 64 -o ${sample_id}.filter_pe.filtered.sorted.bam ${sample_id}.filter_pe.filtered.bam
    samtools index -@ 64 ${sample_id}.filter_pe.filtered.sorted.bam
    rm ${sample_id}.filter_pe.filtered.bam ${sample_id}.filter_pe.unique.sorted.filtered.sam
    samtools idxstats ${sample_id}.filter_pe.filtered.sorted.bam > ${sample_id}.filter_pe.filtered.idxstats.txt

    if [[ -f "${sample_id}.filter_pe.unique.sorted.wrongmap.r1.fastq" ]]; then
        mv ${sample_id}.filter_pe.unique.sorted.wrongmap.r1.fastq ${sample_id}.filter_pe.wrongmap.r1.fastq
    else
        touch ${sample_id}.filter_pe.wrongmap.r1.fastq
    fi
    pigz --best -p 64 ${sample_id}.filter_pe.wrongmap.r1.fastq

    if [[ -f "${sample_id}.filter_pe.unique.sorted.wrongmap.r2.fastq" ]]; then
        mv ${sample_id}.filter_pe.unique.sorted.wrongmap.r2.fastq ${sample_id}.filter_pe.wrongmap.r2.fastq
    else
        touch ${sample_id}.filter_pe.wrongmap.r2.fastq
    fi
    pigz --best -p 64 ${sample_id}.filter_pe.wrongmap.r2.fastq
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

process extract_barcodes_pe {
    label 'process_high'

    input:
    tuple val(sample_id), path(bam), path(bai), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.filter_pe.barcodes.txt"), emit: ch_bwa_pe_barcodes

    script:
    """
    ${projectDir}/scripts/extract_barcodes_pe.R -b ${bam} -a ${barcode} -p ${params.barcode_template} -m ${params.barcode_marker}
    mv ${sample_id}.filter_pe.filtered.sorted.barcodes.txt ${sample_id}.filter_pe.barcodes.txt
    """
}