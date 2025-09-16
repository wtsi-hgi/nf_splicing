workflow create_splicing_counts {
    take:
    ch_sample

    main:
    /* -- 1. classify novel junctions -- */
    ch_novel_junctions = ch_sample.map {sample_id, canonical_barcodes, novel_junctions, exon_pos-> 
                                            tuple(sample_id, novel_junctions, exon_pos)}
    CLASSIFY_NOVEL_JUNCTIONS(ch_novel_junctions)
    ch_classified_junctions = CLASSIFY_NOVEL_JUNCTIONS.out.ch_classified_junctions

    /* -- 2. create splicing counts -- */

}

process CLASSIFY_NOVEL_JUNCTIONS {
    label 'process_medium'

    input:
    tuple val(sample_id), path(novel_junctions), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.txt"), emit: ch_classified_junctions

    script:
    """
    python classify_novel_junctions.py --bed_file ${novel_junctions} \
                                       --exon_pos ${exon_pos} \
                                       --cluster_tol ${params.classify_cluster_tol} \
                                       --junc_cov ${params.classify_min_cov} \
                                       --min_overlap ${params.classify_min_overlap} \
                                       --output_prefix ${sample_id}
    """
}
