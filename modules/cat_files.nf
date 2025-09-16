process cat_canonical_barcodes {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(filter_se), path(filter_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_barcodes.txt"), emit: ch_canonical_barcodes
    
    script:
    """
    head -n 1 ${filter_se} > header.txt
    tail -n +2 ${filter_se} > se.txt
    tail -n +2 ${filter_pe} > pe.txt
    cat header.txt se.txt pe.txt > ${sample_id}.canonical_barcodes.txt
    rm header.txt se.txt pe.txt
    """
}

process cat_novel_barcodes {
    label 'process_single'
    
    input:
    tuple val(sample_id), path(map_se), path(map_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_barcodes.txt"), emit: ch_novel_barcodes
    
    script:
    """
    head -n 1 ${map_se} > header.txt
    tail -n +2 ${map_se} > se.txt
    tail -n +2 ${map_pe} > pe.txt
    cat header.txt se.txt pe.txt > ${sample_id}.novel_barcodes.txt
    rm header.txt se.txt pe.txt
    """
}

process cat_beds {
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