workflow create_splicing_counts {
    take:
    ch_sample

    main:
    /* -- 1. classify novel junctions -- */
    ch_novel_junctions = ch_sample.map {sample_id, barcode, canonical_barcodes, novel_junctions, exon_pos-> 
                                            tuple(sample_id, novel_junctions, exon_pos)}
    CLASSIFY_NOVEL_JUNCTIONS(ch_novel_junctions)
    ch_classified_junctions = CLASSIFY_NOVEL_JUNCTIONS.out.ch_classified_junctions

    /* -- 2. create splicing counts -- */
    ch_input_count = ch_sample.map {sample_id, barcode, canonical_barcodes, novel_junctions, exon_pos-> 
                                        tuple(sample_id, barcode, canonical_barcodes)}
                              .join(ch_classified_junctions)
    CREATE_SPLICING_COUNTS(ch_input_count)
    ch_splicing_counts = CREATE_SPLICING_COUNTS.out.ch_splicing_counts

    emit:
    ch_classified_junctions
    ch_splicing_counts
}

process CLASSIFY_NOVEL_JUNCTIONS {
    label 'process_single_dynamic_mem'

    memory {
        def file_size = novel_junctions.size()
        def mem = file_size <= 1_000_000 ? 2 :
                  file_size <= 10_000_000 ? 8 :
                  file_size <= 100_000_000 ? 16 : 32
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.classified_junctions.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(novel_junctions), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.tsv"), emit: ch_classified_junctions

    script:
    """
    python ${projectDir}/scripts/classify_novel_junctions.py --lib_type ${params.library} \
                                                             --bed_file ${novel_junctions} \
                                                             --exon_pos ${exon_pos} \
                                                             --cluster_tol ${params.classify_cluster_tol} \
                                                             --junc_cov ${params.classify_min_cov} \
                                                             --min_overlap ${params.classify_min_overlap} \
                                                             --output_prefix ${sample_id}
    """
}

process CREATE_SPLICING_COUNTS {
    label 'process_single_dynamic_mem'

    memory {
        def file_size = canonical_barcodes.size()
        def mem = file_size <= 10_000_000 ? 2 :
                  file_size <= 100_000_000 ? 4 :
                  file_size <= 1000_000_000 ? 6 : 8
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/splicing_counts/", pattern: "*.splicing_counts.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), path(canonical_barcodes), path(classified_junctions)

    output:
    tuple val(sample_id), path("${sample_id}.splicing_counts.tsv"), emit: ch_splicing_counts

    script:
    """
    python ${projectDir}/scripts/create_splicing_counts.py --canonical_file ${canonical_barcodes} \
                                                           --novel_file ${classified_junctions} \
                                                           --barcode_file ${barcode} \
                                                           --output_prefix ${sample_id}
    """
}