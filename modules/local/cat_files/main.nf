process CAT_CANONICAL_BARCODES {
    label 'process_single'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", mode: "copy", overwrite: true
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(can_barcode_se), path(can_barcode_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.canonical_barcodes.tsv.gz"), emit: ch_canonical_barcodes
    
    script:
    """
    zcat ${can_barcode_se} | head -n 1 > header.tsv
    zcat ${can_barcode_se} | tail -n +2 > se.tsv
    zcat ${can_barcode_pe} | tail -n +2 > pe.tsv
    cat header.tsv se.tsv pe.tsv > ${sample_id}.canonical_barcodes.tsv
    gzip ${sample_id}.canonical_barcodes.tsv
    rm header.tsv se.tsv pe.tsv
    """
}

process CAT_NOVEL_BARCODES {
    label 'process_single'
    
    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(nov_barcode_se), path(nov_barcode_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.novel_barcodes.tsv.gz"), emit: ch_novel_barcodes
    
    script:
    """
    zcat ${nov_barcode_se} | head -n 1 > header.tsv
    zcat ${nov_barcode_se} | tail -n +2 > se.tsv
    zcat ${nov_barcode_pe} | tail -n +2 > pe.tsv
    cat header.tsv se.tsv pe.tsv > ${sample_id}.novel_barcodes.tsv
    gzip ${sample_id}.novel_barcodes.tsv
    rm header.tsv se.tsv pe.tsv
    """
}

process CAT_BEDS {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bed_se), path(bed_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.junctions.bed"), emit: ch_bed
    
    script:
    """
    cat ${bed_se} ${bed_pe} > ${sample_id}.junctions.bed
    """
}

process CAT_BASE_COVS {
    label 'process_single'

    publishDir "${params.outdir}/splicing_counts/${sample_id}", mode: "copy", overwrite: true
    
    tag "$sample_id"

    input:
    tuple val(sample_id), path(base_cov_se), path(base_cov_pe)
    
    output:
    tuple val(sample_id), path("${sample_id}.base_cov.tsv.gz"), emit: ch_base_cov
    
    script:
    """
    awk 'BEGIN{
            FS="\\t"
            OFS="\\t"
        } NR==FNR {
            a[\$1"_"\$2] = \$3
            next
        }{
            key = \$1"_"\$2
            if (key in a) {
                print \$1, \$2, \$3 + a[key]
            } else {
                print
            }
        }' <(zcat ${base_cov_pe}) <(zcat ${base_cov_se}) | gzip > ${sample_id}.base_cov.tsv.gz
    """    
}
