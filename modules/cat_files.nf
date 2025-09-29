process CAT_CANONICAL_BARCODES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(filter_se), path(filter_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_barcodes.tsv"), emit: ch_canonical_barcodes
    
    script:
    """
    head -n 1 ${filter_se} > header.tsv
    tail -n +2 ${filter_se} > se.tsv
    tail -n +2 ${filter_pe} > pe.tsv
    cat header.tsv se.tsv pe.tsv > ${sample_id}.canonical_barcodes.tsv
    rm header.tsv se.tsv pe.tsv
    """
}

process CAT_NOVEL_BARCODES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(map_se), path(map_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_barcodes.tsv"), emit: ch_novel_barcodes
    
    script:
    """
    head -n 1 ${map_se} > header.tsv
    tail -n +2 ${map_se} > se.tsv
    tail -n +2 ${map_pe} > pe.tsv
    cat header.tsv se.tsv pe.tsv > ${sample_id}.novel_barcodes.tsv
    rm header.tsv se.tsv pe.tsv
    """
}

process CAT_BEDS {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(bed_se), path(bed_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.junctions.bed"), emit: ch_bed
    
    script:
    """
    cat ${bed_se} ${bed_pe} > ${sample_id}.junctions.bed
    """
}