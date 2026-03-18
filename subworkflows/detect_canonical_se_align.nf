include { BWA_ALIGN_SE_READS } from "$projectDir/modules/local/bwa/main"
include { FILTER_SE_READS }    from "$projectDir/modules/local/filter_reads/main"
include { SORT_BWA_SE_BAM }    from "$projectDir/modules/local/sort_bwa_bam/main" 
include { SE_BAM_TO_FASTQ }    from "$projectDir/modules/local/bam_to_fastq/main"

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

    SORT_BWA_SE_BAM(ch_se_filtered_bam)
    ch_se_sorted_bam = SORT_BWA_SE_BAM.out.ch_se_sorted_bam
    ch_se_canonical_stats = SORT_BWA_SE_BAM.out.ch_se_canonical_stats

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
