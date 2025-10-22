workflow detect_canonical_se_align {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos -> 
                                        tuple(sample_id, extended_frags, exon_fasta) }
    BWA_ALIGN_SE_READS(ch_align_input)
    ch_se_unique_bam = BWA_ALIGN_SE_READS.out.ch_se_unique_bam
    ch_se_unmapped = BWA_ALIGN_SE_READS.out.ch_se_unmapped

    /* -- 2. filter reads -- */
    ch_filter_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, exon_pos) }
                               .join(ch_se_unique_bam)

    FILTER_SE_READS(ch_filter_input)
    ch_se_filtered_bam = FILTER_SE_READS.out.ch_se_filtered_bam
    ch_se_wrongmap_bam = FILTER_SE_READS.out.ch_se_wrongmap_bam
    ch_se_canonical_barcodes = FILTER_SE_READS.out.ch_se_canonical_barcodes

    SORT_SE_BAM(ch_se_filtered_bam)
    ch_se_sorted_bam = SORT_SE_BAM.out.ch_se_sorted_bam
    ch_se_canonical_stats = SORT_SE_BAM.out.ch_se_canonical_stats

    SE_BAM_TO_FASTQ(ch_se_wrongmap_bam)
    ch_se_wrongmap = SE_BAM_TO_FASTQ.out.ch_se_wrongmap

    /* -- 3. cat reads -- */
    ch_se_joined = ch_se_unmapped.join(ch_se_wrongmap)
    MERGE_SE_FASTQS(ch_se_joined)
    ch_se_canonical_fail = MERGE_SE_FASTQS.out.ch_se_canonical_fail

    emit:
    ch_se_sorted_bam
    ch_se_canonical_fail
    ch_se_canonical_barcodes
    ch_se_canonical_stats
}

process BWA_ALIGN_SE_READS {
    label 'process_medium'

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

process FILTER_SE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(exon_pos), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.filtered.bam"), emit: ch_se_filtered_bam
    tuple val(sample_id), path("${sample_id}.bwa_se.wrongmap.bam"), emit: ch_se_wrongmap_bam
    tuple val(sample_id), path("${sample_id}.bwa_se.canonical_barcodes.tsv"), emit: ch_se_canonical_barcodes

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
                                                          --output_prefix ${sample_id}.bwa_se \
                                                          --chunk_size 100000 \
                                                          --threads 40
    """
}

process SORT_SE_BAM {
    label 'process_high_memory'

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
    """
}

process SE_BAM_TO_FASTQ {
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.wrongmap.fastq.gz"), emit: ch_se_wrongmap

    script:
    """
    samtools fastq -@ 32 -c 9 -0 ${sample_id}.bwa_se.wrongmap.fastq.gz ${bam}
    rm ${bam}
    """
}

process MERGE_SE_FASTQS {
    label 'process_single'

    input:
    tuple val(sample_id), path(unmapped_fastq), path(wrongmap_fastq)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.canonical_fail.fastq.gz"), emit: ch_se_canonical_fail

    script:
    """
    cat ${unmapped_fastq} ${wrongmap_fastq} > ${sample_id}.bwa_se.canonical_fail.fastq.gz
    """
}
