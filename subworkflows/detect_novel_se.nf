include { HISAT2_ALIGN_SE_READS } from "$projectDir/modules/local/hisat2/main"
include { FIX_SE_READS }          from "$projectDir/modules/local/fix_reads/main"
include { SORT_SE_BAM }           from "$projectDir/modules/local/sort_hisat2_bam/main"
include { EXTRACT_SE_JUNCTIONS }  from "$projectDir/modules/local/regtools/main"

workflow detect_novel_se {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, se_reads, ref_fasta -> 
                                        tuple(sample_id, se_reads, ref_fasta) }
    HISAT2_ALIGN_SE_READS(ch_align_input)
    ch_se_unique_bam = HISAT2_ALIGN_SE_READS.out.ch_se_unique_bam
    ch_se_novel_stats = HISAT2_ALIGN_SE_READS.out.ch_se_novel_stats

    /* -- 2. fix alignments -- */
    ch_fix_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, se_reads, ref_fasta -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, ref_fasta) }
                            .join(ch_se_unique_bam)
    
    FIX_SE_READS(ch_fix_input)
    ch_se_fixed_bam = FIX_SE_READS.out.ch_se_fixed_bam
    ch_se_novel_barcodes = FIX_SE_READS.out.ch_se_novel_barcodes

    SORT_SE_BAM(ch_se_fixed_bam)
    ch_se_sorted_bam = SORT_SE_BAM.out.ch_se_sorted_bam

    /* -- 3. extract junctions -- */
    EXTRACT_SE_JUNCTIONS(ch_se_sorted_bam)
    ch_se_junctions = EXTRACT_SE_JUNCTIONS.out.ch_se_junctions

    emit:
    ch_se_sorted_bam
    ch_se_novel_barcodes
    ch_se_junctions
    ch_se_novel_stats
}
