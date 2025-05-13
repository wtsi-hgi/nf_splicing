process publish_canonical_barcodes {
    label 'process_single'
    
    publishDir "${params.outdir}/extracted_barcodes/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(output_file)
    
    output:
    tuple val(sample_id), path(output_file), emit: ch_publish
    
    script:
    """
    filename=\$(basename ${output_file})
    echo "Publish \$filename for sample: ${sample_id}"
    """
}

process publish_novel_barcodes {
    label 'process_single'
    
    publishDir "${params.outdir}/extracted_barcodes/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(output_file)
    
    output:
    tuple val(sample_id), path(output_file), emit: ch_publish
    
    script:
    """
    filename=\$(basename ${output_file})
    echo "Publish \$filename for sample: ${sample_id}"
    """
}