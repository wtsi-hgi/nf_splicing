workflow detect_canonical_pe {
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
    ch_pe_wrongmap = FILTER_PE_READS.out.ch_pe_wrongmap
    ch_pe_canonical_barcodes = FILTER_PE_READS.out.ch_pe_canonical_barcodes
    ch_pe_canonical_stats = FILTER_PE_READS.out.ch_pe_canonical_stats

    /* -- 3. cat reads -- */
    ch_pe_joined = ch_pe_unmapped.join(ch_pe_wrongmap)
    MERGE_PE_FASTQS(ch_pe_joined)
    ch_pe_canonical_fail = MERGE_PE_FASTQS.out.ch_pe_canonical_fail

    emit:
    ch_pe_filtered_bam
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

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.barcodes.tsv", mode: "copy", overwrite: true

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
    samtools sort -@ 40 -o ${sample_id}.bwa_pe.filtered.sorted.bam ${sample_id}.bwa_pe.filtered.bam
    samtools index -@ 40 ${sample_id}.bwa_pe.filtered.sorted.bam
    rm ${sample_id}.bwa_pe.filtered.bam
    samtools idxstats ${sample_id}.bwa_pe.filtered.sorted.bam > ${sample_id}.bwa_pe.canonical_stats.tsv

    samtools fastq -@ 40 -c 9 -1 ${sample_id}.bwa_pe.wrongmap.r1.fastq.gz -2 ${sample_id}.bwa_pe.wrongmap.r2.fastq.gz -n ${sample_id}.bwa_pe.wrongmap.bam

    rm ${sample_id}.bwa_pe.filtered.bam
    rm ${sample_id}.bwa_pe.wrongmap.bam
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
