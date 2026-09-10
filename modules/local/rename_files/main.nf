process RENAME_CANONICAL_BARCODES {
    label 'process_single'
    
    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv.gz", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(filter_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_barcodes.tsv.gz"), emit: ch_canonical_barcodes
    
    script:
    """
    mv ${filter_se} ${sample_id}.canonical_barcodes.tsv.gz
    """
}

process RENAME_NOVEL_BARCODES {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.novel_barcodes.tsv.gz", mode: "copy", overwrite: true
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(map_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_barcodes.tsv.gz"), emit: ch_novel_barcodes
    
    script:
    """
    mv ${map_se} ${sample_id}.novel_barcodes.tsv.gz
    """
}
