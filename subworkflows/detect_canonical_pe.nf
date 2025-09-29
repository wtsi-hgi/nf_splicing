workflow detect_canonical_pe {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                        tuple(sample_id, not_combined_1, not_combined_2, exon_fasta) }
    BWA_ALIGN_PE_READS(ch_align_input)
    ch_bwa_pe_unmapped = BWA_ALIGN_PE_READS.out.ch_bwa_pe_unmapped
    ch_bwa_pe_bam = BWA_ALIGN_PE_READS.out.ch_bwa_pe_bam

    /* -- 2. filter reads -- */
    ch_filter_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, exon_pos) }
                               .join(ch_bwa_pe_bam)    
    FILTER_PE_READS(ch_filter_input)
    ch_bwa_pe_wrongmap = FILTER_PE_READS.out.ch_bwa_pe_wrongmap
    ch_bwa_pe_filtered = FILTER_PE_READS.out.ch_bwa_pe_filtered
    ch_bwa_pe_filtered_idxstats = FILTER_PE_READS.out.ch_bwa_pe_filtered_idxstats
    ch_bwa_pe_barcodes = FILTER_PE_READS.out.ch_bwa_pe_barcodes

    /* -- 3. cat reads -- */
    ch_bwa_pe_joined = ch_bwa_pe_unmapped.join(ch_bwa_pe_wrongmap)
    MERGE_PE_FASTQS(ch_bwa_pe_joined)
    ch_bwa_pe_fail = MERGE_PE_FASTQS.out.ch_bwa_pe_fail

    emit:
    ch_bwa_pe_filtered
    ch_bwa_pe_filtered_idxstats
    ch_bwa_pe_fail
    ch_bwa_pe_barcodes
}

process BWA_ALIGN_PE_READS {
    label 'process_medium'

    input:
    tuple val(sample_id), path(not_combined_1), path(not_combined_2), path(exon_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.filter_pe.unmapped.r1.fastq.gz"), path("${sample_id}.filter_pe.unmapped.r2.fastq.gz"), emit: ch_bwa_pe_unmapped
    tuple val(sample_id), path("${sample_id}.filter_pe.unique.bam"), emit: ch_bwa_pe_bam

    script:
    """
    bwa index ${exon_fasta}
    bwa mem -t 32 -O ${params.bwa_gap_open} \
                  -E ${params.bwa_gap_ext} \
                  -L ${params.bwa_clip} \
                  ${exon_fasta} ${not_combined_1} ${not_combined_2} > ${sample_id}.filter_pe.sam
    samtools fastq -@ 32 -F 2 -c 9 -1 ${sample_id}.filter_pe.unmapped.r1.fastq.gz -2 ${sample_id}.filter_pe.unmapped.r2.fastq.gz -n ${sample_id}.filter_pe.sam
    samtools view -b -f 2 -F 256 -F 2048 ${sample_id}.filter_pe.sam > ${sample_id}.filter_pe.unique.bam
    rm ${sample_id}.filter_pe.sam
    """
}

process FILTER_PE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.barcodes.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(exon_pos), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.filter_pe.wrongmap.r1.fastq.gz"), path("${sample_id}.filter_pe.wrongmap.r2.fastq.gz"), emit: ch_bwa_pe_wrongmap
    tuple val(sample_id), path("${sample_id}.filter_pe.filtered.sorted.bam"), path("${sample_id}.filter_pe.filtered.sorted.bam.bai"), emit: ch_bwa_pe_filtered
    tuple val(sample_id), path("${sample_id}.filter_pe.filtered.idxstats.tsv"), emit: ch_bwa_pe_filtered_idxstats
    tuple val(sample_id), path("${sample_id}.filter_pe.barcodes.tsv"), emit: ch_bwa_pe_barcodes

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
                                                          --output_prefix ${sample_id}.filter_pe \
                                                          --chunk_size 100000 \
                                                          --threads 40
    samtools sort -@ 40 -o ${sample_id}.filter_pe.filtered.sorted.bam ${sample_id}.filter_pe.filtered.bam
    samtools index -@ 40 ${sample_id}.filter_pe.filtered.sorted.bam
    rm ${sample_id}.filter_pe.filtered.bam
    samtools idxstats ${sample_id}.filter_pe.filtered.sorted.bam > ${sample_id}.filter_pe.filtered.idxstats.tsv

    samtools fastq -@ 40 -c 9 -1 ${sample_id}.filter_pe.wrongmap.r1.fastq.gz -2 ${sample_id}.filter_pe.wrongmap.r2.fastq.gz -n ${sample_id}.filter_pe.wrongmap.bam
    """
}

process MERGE_PE_FASTQS {
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
