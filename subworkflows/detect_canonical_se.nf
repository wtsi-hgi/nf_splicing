workflow detect_canonical_se {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos -> 
                                        tuple(sample_id, extended_frags, exon_fasta) }
    bwa_align_se_reads(ch_align_input)
    ch_bwa_se_unmapped = bwa_align_se_reads.out.ch_bwa_se_unmapped
    ch_bwa_se_bam = bwa_align_se_reads.out.ch_bwa_se_bam

    /* -- 2. filter reads -- */
    ch_filter_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, exon_pos) }
                               .join(ch_bwa_se_bam)
    filter_se_reads(ch_filter_input)
    ch_bwa_se_wrongmap = filter_se_reads.out.ch_bwa_se_wrongmap
    ch_bwa_se_filtered = filter_se_reads.out.ch_bwa_se_filtered
    ch_bwa_se_filtered_idxstats = filter_se_reads.out.ch_bwa_se_filtered_idxstats
    ch_bwa_se_barcodes = filter_se_reads.out.ch_bwa_se_barcodes

    /* -- 3. cat reads -- */
    ch_bwa_se_joined = ch_bwa_se_unmapped.join(ch_bwa_se_wrongmap)
    merge_se_fastqs(ch_bwa_se_joined)
    ch_bwa_se_fail = merge_se_fastqs.out.ch_bwa_se_fail

    emit:
    ch_bwa_se_filtered
    ch_bwa_se_filtered_idxstats
    ch_bwa_se_fail
    ch_bwa_se_barcodes
}

process bwa_align_se_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(extended_frags), path(exon_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.unmapped.fastq.gz"), emit: ch_bwa_se_unmapped
    tuple val(sample_id), path("${sample_id}.filter_se.unique.bam"), emit: ch_bwa_se_bam

    script:
    """
    bwa index ${exon_fasta}
    bwa mem -t 32 -O ${params.bwa_gap_open} \
                  -E ${params.bwa_gap_ext} \
                  -L ${params.bwa_clip} \
                  ${exon_fasta} ${extended_frags} > ${sample_id}.filter_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.filter_se.unmapped.fastq.gz ${sample_id}.filter_se.sam
    samtools view -@ 32 -b -F 4 -F 256 -F 2048 ${sample_id}.filter_se.sam > ${sample_id}.filter_se.unique.bam
    rm ${sample_id}.filter_se.sam
    """
}

process filter_se_reads {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.barcodes.txt", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(exon_pos), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.filter_se.wrongmap.fastq.gz"), emit: ch_bwa_se_wrongmap
    tuple val(sample_id), path("${sample_id}.filter_se.filtered.sorted.bam"), path("${sample_id}.filter_se.filtered.sorted.bam.bai"), emit: ch_bwa_se_filtered
    tuple val(sample_id), path("${sample_id}.filter_se.filtered.idxstats.txt"), emit: ch_bwa_se_filtered_idxstats
    tuple val(sample_id), path("${sample_id}.filter_se.barcodes.txt"), emit: ch_bwa_se_barcodes

    script:
    """
    python ${projectDir}/scripts/process_canonical_bam.py --bam_file ${bam} \
                                                          --barcode_file ${barcode} \
                                                          --exon_pos ${exon_pos} \
                                                          --read_type se \
                                                          --soft_clip ${params.filter_softclip_base} \
                                                          --barcode_up ${barcode_up} \
                                                          --barcode_down ${barcode_down} \
                                                          --barcode_check \
                                                          --barcode_temp ${barcode_temp} \
                                                          --output_prefix ${sample_id}.filter_se \
                                                          --chunk_size 100000 \
                                                          --threads 40
    samtools sort -@ 40 -o ${sample_id}.filter_se.filtered.sorted.bam ${sample_id}.filter_se.filtered.bam
    samtools index -@ 40 ${sample_id}.filter_se.filtered.sorted.bam
    rm ${sample_id}.filter_se.filtered.bam
    samtools idxstats ${sample_id}.filter_se.filtered.sorted.bam > ${sample_id}.filter_se.filtered.idxstats.txt

    samtools fastq -@ 40 -c 9 -0 ${sample_id}.filter_se.wrongmap.fastq.gz ${sample_id}.filter_se.wrongmap.bam
    """
}

process merge_se_fastqs {
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
