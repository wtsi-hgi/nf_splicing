process RENAME_CANONICAL_BARCODES {
    label 'process_single'
    
    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", mode: "copy", overwrite: true

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

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true
    
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

process RENAME_BEDS {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bed_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.junctions.bed"), emit: ch_junctions
    
    script:
    """
    mv ${bed_se} ${sample_id}.junctions.bed
    """
}

process RENAME_BASE_COVS {
    label 'process_single'

    publishDir "${params.outdir}/splicing_counts/${sample_id}", mode: "copy", overwrite: true
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(base_cov_se)
    
    output:
    tuple val(sample_id), path("${sample_id}.base_cov.tsv.gz"), emit: ch_base_cov
    
    script:
    """
    mv ${base_cov_se} ${sample_id}.base_cov.tsv.gz
    """
}
