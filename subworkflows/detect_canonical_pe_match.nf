include { MATCH_PE_READS } from "$projectDir/modules/local/match_reads/main"

workflow detect_canonical_pe_match {
    take:
    ch_sample

    main:
    /* -- 1. match reads -- */
    MATCH_PE_READS(ch_sample)
    ch_pe_canonical_fail = MATCH_PE_READS.out.ch_pe_canonical_fail
    ch_pe_canonical_barcodes = MATCH_PE_READS.out.ch_pe_canonical_barcodes
    ch_pe_unknown = MATCH_SE_READS.out.ch_pe_unknown

    /* -- 2. get stats -- */
    GET_STATS(ch_pe_canonical_barcodes)
    ch_pe_canonical_stats = GET_STATS.out.ch_pe_canonical_stats

    emit:
    ch_pe_unknown
    ch_pe_canonical_fail
    ch_pe_canonical_barcodes
    ch_pe_canonical_stats
}

process GET_STATS {
    label 'process_single'

    input:
    tuple val(sample_id), path(pe_canonical_barcodes)

    output:
    tuple val(sample_id), path("${sample_id}.match_pe.canonical_stats.tsv"), emit: ch_pe_canonical_stats

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
    }' ${pe_canonical_barcodes} > ${sample_id}.match_pe.canonical_stats.tsv
    """
}
