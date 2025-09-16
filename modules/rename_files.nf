process RENAME_CANONICAL_BARCODES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(filter_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_barcodes.txt"), emit: ch_canonical_barcodes
    
    script:
    """
    mv ${filter_se} ${sample_id}.canonical_barcodes.txt
    """
}

process RENAME_NOVEL_BARCODES {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(map_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_barcodes.txt"), emit: ch_novel_barcodes
    
    script:
    """
    mv ${map_se} ${sample_id}.novel_barcodes.txt
    """
}
