include { CLASSIFY_NOVEL_JUNCTIONS } from "$projectDir/modules/local/classify_junctions/main"
include { CREATE_SPLICING_COUNTS }   from "$projectDir/modules/local/create_counts/main"
include { CALCULATE_SSU }            from "$projectDir/modules/local/calculate_ssu/main"

workflow create_splicing_counts {
    take:
    ch_sample

    main:
    /* -- 1. classify novel junctions -- */
    ch_novel_junctions = ch_sample.map { sample_id, barcode, canonical_barcodes, novel_junctions, base_cov, exon_pos-> 
                                            tuple(sample_id, novel_junctions, exon_pos) }
    CLASSIFY_NOVEL_JUNCTIONS(ch_novel_junctions)
    ch_classified_junctions = CLASSIFY_NOVEL_JUNCTIONS.out.ch_classified_junctions

    /* -- 2. create splicing counts -- */
    ch_input_count = ch_sample.map { sample_id, barcode, canonical_barcodes, novel_junctions, base_cov, exon_pos-> 
                                        tuple(sample_id, barcode, canonical_barcodes) }
                              .join(ch_classified_junctions)
    CREATE_SPLICING_COUNTS(ch_input_count)
    ch_splicing_counts = CREATE_SPLICING_COUNTS.out.ch_splicing_counts

    /* -- 3. create splicing site usage (SSU) -- */
    ch_input_count = ch_sample.map { sample_id, barcode, canonical_barcodes, novel_junctions, base_cov, exon_pos-> 
                                        tuple(sample_id, base_cov, exon_pos) }
                              .join(ch_splicing_counts)
    CALCULATE_SSU(ch_input_count)
    ch_splicing_ssu = CALCULATE_SSU.out.ch_splicing_ssu

    emit:
    ch_classified_junctions
    ch_splicing_counts
    ch_splicing_ssu
}

