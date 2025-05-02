process idxstats_get_values {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(idxstats_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.idxstats.txt"), emit: ch_idxstats
    
    script:
    """
    awk -F'\\t' -v OFS='\\t' '{
        if (\$1!="*") {print \$1, \$3}
    }' ${idxstats_se} | sort -k1,1 >  ${sample_id}.idxstats.txt
    """
}

process idxstats_add_values {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(idxstats_se), path(idxstats_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.idxstats.txt"), emit: ch_idxstats
    
    script:
    """
    cat ${idxstats_se} ${idxstats_pe} | awk -F'\\t' -v OFS='\\t' '{
        if (\$1!="*") {
            if (a[\$1]) {
                a[\$1] = a[\$1] + \$3
            } else {
                a[\$1] = \$3
            }
        }
    }END{
        for (i in a) {
            print i, a[i]
        }
    }' | sort -k1,1 > ${sample_id}.idxstats.txt
    """
}
