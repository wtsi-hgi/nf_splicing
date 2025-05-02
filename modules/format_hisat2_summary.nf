process hisat2_summary_get_values {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(hisat2_summary_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.hisat2_summary.txt"), emit: ch_hisat2_summary
    
    script:
    """
    awk -F' ' -v OFS='\\t' '
        /Total reads:/ {total = \$NF}
        /Aligned 1 time:/ {map = \$4}
        END {print "Total", total
             print "Mapped", map}' ${hisat2_summary_se} | sort -k1,1 > ${sample_id}.hisat2_summary.txt
    """
}

process hisat2_summary_add_values {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(hisat2_summary_se), path(hisat2_summary_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.hisat2_summary.txt"), emit: ch_hisat2_summary
    
    script:
    """
    awk -F' ' -v OFS='\\t' '
        /Total reads:/ {total = \$NF}
        /Aligned 1 time:/ {map = \$4}
        END {print "Total", total
             print "Mapped", map}' ${hisat2_summary_se} | sort -k1,1 > se.txt
    awk -F' ' -v OFS='\\t' '
        /Total pairs:/ {total = \$NF}
        /Aligned concordantly 1 time:/ {map = \$5}
        END {print "Total", total
             print "Mapped", map}' ${hisat2_summary_pe} | sort -k1,1 > pe.txt
    awk -F'\\t' -v OFS='\\t' '
        NR==FNR {a[\$1] = \$2; next}
        {print \$1, a[\$1] + \$2}' se.txt pe.txt > ${sample_id}.hisat2_summary.txt
    rm se.txt pe.txt
    """
}
