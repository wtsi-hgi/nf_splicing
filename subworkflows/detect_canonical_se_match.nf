workflow detect_canonical_se_match {
    take:
    ch_sample

    main:
    /* -- 1. match reads -- */
    MATCH_SE_READS(ch_sample)
    ch_se_canonical_fail = MATCH_SE_READS.out.ch_se_canonical_fail
    ch_se_canonical_barcodes = MATCH_SE_READS.out.ch_se_canonical_barcodes
    ch_se_unknown = MATCH_SE_READS.out.ch_se_unknown

    /* -- 2. get stats -- */
    GET_STATS(ch_se_canonical_barcodes)
    ch_se_canonical_stats = GET_STATS.out.ch_se_canonical_stats

    emit:
    ch_se_canonical_fail
    ch_se_canonical_barcodes
    ch_se_unknown
    ch_se_canonical_stats
}

process MATCH_SE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(extended_frags), path(ch_hisat2_ref)

    output:
    tuple val(sample_id), path("${sample_id}.match_se.canonical_fail.fastq.gz"), emit: ch_se_canonical_fail
    tuple val(sample_id), path("${sample_id}.match_se.canonical_barcodes.tsv"), emit: ch_se_canonical_barcodes
    tuple val(sample_id), path("${sample_id}.match_se.unknown.tsv"), emit: ch_se_unknown

    script:
    """
    python ${projectDir}/scripts/process_canonical_fastq.py --reads ${extended_frags} \
                                                            --read_type se \
                                                            --ref_file ${ch_hisat2_ref} \
                                                            --barcode_file ${barcode} \
                                                            --barcode_up ${barcode_up} \
                                                            --barcode_down ${barcode_down} \
                                                            --barcode_check \
                                                            --barcode_temp ${barcode_temp} \
                                                            --output_prefix ${sample_id}.match_se \
                                                            --chunk_size 100000 \
                                                            --threads 40
    """
}

process GET_STATS {
    label 'process_single'

    input:
    tuple val(sample_id), path(se_canonical_barcodes)

    output:
    tuple val(sample_id), path("${sample_id}.match_se.canonical_stats.tsv"), emit: ch_se_canonical_stats

    script:
    """
    awk -F'\\t' -v OFS='\\t' 'BEGIN{cin = 0; csk = 0}{
        if(\$5 == "exon_inclusion") {
            cin = cin + \$4
        } else if(\$5 == "exon_skipping") {
            csk = csk + \$4
        }
    }END{
        print "exon_inclusion", cin;
        print "exon_skipping", csk
    }' ${se_canonical_barcodes} > ${sample_id}.match_se.canonical_stats.tsv
    """
}
