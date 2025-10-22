workflow detect_canonical_pe_align {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                        tuple(sample_id, not_combined_1, not_combined_2, exon_fasta) }
    BWA_ALIGN_PE_READS(ch_align_input)
    ch_pe_unique_bam = BWA_ALIGN_PE_READS.out.ch_pe_unique_bam
    ch_pe_unmapped = BWA_ALIGN_PE_READS.out.ch_pe_unmapped

    /* -- 2. filter reads -- */
    ch_filter_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, exon_pos) }
                               .join(ch_pe_unique_bam)    
    FILTER_PE_READS(ch_filter_input)
    ch_pe_filtered_bam = FILTER_PE_READS.out.ch_pe_filtered_bam
    ch_pe_wrongmap_bam = FILTER_PE_READS.out.ch_pe_wrongmap_bam
    ch_pe_canonical_barcodes = FILTER_PE_READS.out.ch_pe_canonical_barcodes

    SORT_PE_BAM(ch_pe_filtered_bam)
    ch_pe_sorted_bam = SORT_PE_BAM.out.ch_pe_sorted_bam
    ch_pe_canonical_stats = SORT_PE_BAM.out.ch_pe_canonical_stats

    PE_BAM_TO_FASTQ(ch_pe_wrongmap_bam)
    ch_pe_wrongmap = PE_BAM_TO_FASTQ.out.ch_pe_wrongmap

    /* -- 3. cat reads -- */
    ch_pe_joined = ch_pe_unmapped.join(ch_pe_wrongmap)
    MERGE_PE_FASTQS(ch_pe_joined)
    ch_pe_canonical_fail = MERGE_PE_FASTQS.out.ch_pe_canonical_fail

    emit:
    ch_pe_sorted_bam
    ch_pe_canonical_fail
    ch_pe_canonical_barcodes
    ch_pe_canonical_stats
}

process BWA_ALIGN_PE_READS {
    label 'process_medium'

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

process FILTER_PE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(exon_pos), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.filtered.sorted.bam"), path("${sample_id}.bwa_pe.filtered.sorted.bam.bai"), emit: ch_pe_filtered_bam
    tuple val(sample_id), path("${sample_id}.bwa_pe.wrongmap.r1.fastq.gz"), path("${sample_id}.bwa_pe.wrongmap.r2.fastq.gz"), emit: ch_pe_wrongmap
    tuple val(sample_id), path("${sample_id}.bwa_pe.canonical_barcodes.tsv"), emit: ch_pe_canonical_barcodes
    tuple val(sample_id), path("${sample_id}.bwa_pe.canonical_stats.tsv"), emit: ch_pe_canonical_stats


    script:
    """
    python ${projectDir}/scripts/process_canonical_bam.py --bam_file ${bam} \
                                                          --barcode_file ${barcode} \
                                                          --exon_pos ${exon_pos} \
                                                          --read_type pe \
                                                          --soft_clip ${params.filter_softclip_base} \
                                                          --barcode_up ${barcode_up} \
                                                          --barcode_down ${barcode_down} \
                                                          --barcode_check \
                                                          --barcode_temp ${barcode_temp} \
                                                          --output_prefix ${sample_id}.bwa_pe \
                                                          --chunk_size 100000 \
                                                          --threads 40
    """
}

process SORT_PE_BAM {
    label 'process_high_memory'

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
    """
}

process PE_BAM_TO_FASTQ {
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.wrongmap.r1.fastq.gz"), path("${sample_id}.bwa_pe.wrongmap.r2.fastq.gz"), emit: ch_pe_wrongmap

    script:
    """
    samtools fastq -@ 32 -c 9 -1 ${sample_id}.bwa_pe.wrongmap.r1.fastq.gz -2 ${sample_id}.bwa_pe.wrongmap.r2.fastq.gz -n ${bam}
    rm ${bam}
    """
}

process MERGE_PE_FASTQS {
    label 'process_single'

    input:
    tuple val(sample_id), path(unmapped_fastq_r1), path(unmapped_fastq_r2), path(wrongmap_fastq_r1), path(wrongmap_fastq_r2)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.canonical_fail.r1.fastq.gz"), path("${sample_id}.bwa_pe.canonical_fail.r2.fastq.gz"), emit: ch_pe_canonical_fail

    script:
    """
    cat ${unmapped_fastq_r1} ${wrongmap_fastq_r1} > ${sample_id}.bwa_pe.canonical_fail.r1.fastq.gz
    cat ${unmapped_fastq_r2} ${wrongmap_fastq_r2} > ${sample_id}.bwa_pe.canonical_fail.r2.fastq.gz
    """
}
