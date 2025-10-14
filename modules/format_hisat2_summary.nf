process HISAT2_SUMMARY_GET_VALUES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(hisat2_summary_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_stats.tsv"), emit: ch_novel_stats
    
    script:
    """
    awk -F' ' -v OFS='\\t' '
        /Total reads:/ {total = \$NF}
        /Aligned 1 time:/ {map = \$4}
        END {print "Total", total
             print "Mapped", map}' ${hisat2_summary_se} | sort -k1,1 > ${sample_id}.novel_stats.tsv
    """
}

process HISAT2_SUMMARY_ADD_VALUES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(hisat2_summary_se), path(hisat2_summary_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_stats.tsv"), emit: ch_novel_stats
    
    script:
    """
    awk -F' ' -v OFS='\\t' '
        /Total reads:/ {total = \$NF}
        /Aligned 1 time:/ {map = \$4}
        END {print "Total", total
             print "Mapped", map}' ${hisat2_summary_se} | sort -k1,1 > se.tsv
    awk -F' ' -v OFS='\\t' '
        /Total pairs:/ {total = \$NF}
        /Aligned concordantly 1 time:/ {map = \$5}
        END {print "Total", total
             print "Mapped", map}' ${hisat2_summary_pe} | sort -k1,1 > pe.tsv
    awk -F'\\t' -v OFS='\\t' '
        NR==FNR {a[\$1] = \$2; next}
        {print \$1, a[\$1] + \$2}' se.tsv pe.tsv > ${sample_id}.novel_stats.tsv
    rm se.tsv pe.tsv
    """
}
