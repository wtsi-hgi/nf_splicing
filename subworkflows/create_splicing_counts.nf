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

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.classified_junctions.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(novel_junctions), path(exon_pos)

    memory {
        def file_size = novel_junctions.size()     
        file_size <= 1000000 ? check_max( 2.GB * task.attempt, 'memory') :
        file_size <= 10000000 ? check_max( 8.GB * task.attempt, 'memory') :
        file_size <= 100000000 ? check_max( 16.GB * task.attempt, 'memory') : 
        check_max( 32.GB * task.attempt, 'memory')
    }

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.tsv"), emit: ch_classified_junctions

    script:
    """
    python ${projectDir}/scripts/classify_novel_junctions.py --bed_file ${novel_junctions} \
                                                             --exon_pos ${exon_pos} \
                                                             --cluster_tol ${params.classify_cluster_tol} \
                                                             --junc_cov ${params.classify_min_cov} \
                                                             --min_overlap ${params.classify_min_overlap} \
                                                             --output_prefix ${sample_id}
    """
}

process CREATE_SPLICING_COUNTS {
    label 'process_single'

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