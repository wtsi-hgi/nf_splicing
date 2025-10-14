process STATS_GET_VALUES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(se_stats)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_stats.tsv"), emit: ch_canonical_stats
    
    script:
    if (params.canonical_method == 'match') {
        """
        mv ${se_stats} ${sample_id}.canonical_stats.tsv
        """
    } else {
        """
        awk -F'\\t' -v OFS='\\t' '{
            if (\$1!="*") {print \$1, \$3}
        }' ${se_stats} | sort -k1,1 >  ${sample_id}.canonical_stats.tsv
        
        rm ${se_stats}
        """
    }
}

process STATS_ADD_VALUES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(se_stats), path(pe_stats)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_stats.tsv"), emit: ch_canonical_stats
    
    script:
    if (params.canonical_method == 'match') {
        """
        cat ${se_stats} ${pe_stats} | awk -F'\\t' -v OFS='\\t' '{
            if (a[\$1]) {
                a[\$1] = a[\$1] + \$2
            } else {
                a[\$1] = \$2
            }
        }END{
            for (i in a) {
                print i, a[i]
            }
        }' | sort -k1,1 > ${sample_id}.canonical_stats.tsv

        rm ${se_stats} ${pe_stats}
        """    
    } else {
        """
        cat ${se_stats} ${pe_stats} | awk -F'\\t' -v OFS='\\t' '{
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
        }' | sort -k1,1 > ${sample_id}.canonical_stats.tsv

        rm ${se_stats} ${pe_stats}
        """
    }
}
