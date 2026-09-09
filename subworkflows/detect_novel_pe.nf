include { HISAT2_ALIGN_PE_READS } from "$projectDir/modules/local/hisat2/main"
include { FIX_PE_READS }          from "$projectDir/modules/local/fix_reads/main"
include { SORT_PE_BAM }           from "$projectDir/modules/local/sort_hisat2_bam/main"
include { EXTRACT_PE_JUNCTIONS }  from "$projectDir/modules/local/regtools/main"
include { BASE_COV }              from "$projectDir/modules/local/base_cov/main"

workflow detect_novel_pe {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, read1, read2, ref_fasta -> 
                                        tuple(sample_id, read1, read2, ref_fasta) }
    HISAT2_ALIGN_PE_READS(ch_align_input)
    ch_pe_unique_bam = HISAT2_ALIGN_PE_READS.out.ch_pe_unique_bam                               
    ch_pe_novel_stats = HISAT2_ALIGN_PE_READS.out.ch_pe_novel_stats


    /* -- 2. fix alignments -- */
    ch_fix_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, read1, read2, ref_fasta -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, ref_fasta) }
                            .join(ch_pe_unique_bam)

    FIX_PE_READS(ch_fix_input)
    ch_pe_fixed_bam = FIX_PE_READS.out.ch_pe_fixed_bam
    ch_pe_novel_barcodes = FIX_PE_READS.out.ch_pe_novel_barcodes

    SORT_PE_BAM(ch_pe_fixed_bam)
    ch_pe_sorted_bam = SORT_PE_BAM.out.ch_pe_sorted_bam
    
    /* -- 3. extract junctions -- */
    EXTRACT_PE_JUNCTIONS(ch_pe_sorted_bam)
    ch_pe_junctions = EXTRACT_PE_JUNCTIONS.out.ch_pe_junctions

    /* -- 4. get base coverage -- */
    BASE_COV(ch_pe_sorted_bam)
    ch_pe_base_cov = BASE_COV.out.ch_base_cov

    emit:
    ch_pe_sorted_bam
    ch_pe_novel_stats
    ch_pe_novel_barcodes
    ch_pe_junctions
    ch_pe_base_cov
}
