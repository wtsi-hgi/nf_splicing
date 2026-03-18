include { BWA_ALIGN_PE_READS } from "$projectDir/modules/local/bwa/main"
include { FILTER_PE_READS }    from "$projectDir/modules/local/filter_reads/main"
include { SORT_BWA_PE_BAM }    from "$projectDir/modules/local/sort_bwa_bam/main"
include { PE_BAM_TO_FASTQ }    from "$projectDir/modules/local/bam_to_fastq/main"

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

    SORT_BWA_PE_BAM(ch_pe_filtered_bam)
    ch_pe_sorted_bam = SORT_BWA_PE_BAM.out.ch_pe_sorted_bam
    ch_pe_canonical_stats = SORT_BWA_PE_BAM.out.ch_pe_canonical_stats

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
