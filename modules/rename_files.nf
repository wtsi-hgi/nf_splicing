process rename_canonical_barcodes {
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

process rename_novel_barcodes {
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
