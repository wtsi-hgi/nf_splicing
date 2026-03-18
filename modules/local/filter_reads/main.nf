process FILTER_SE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(exon_pos), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_se.filtered.bam"), emit: ch_se_filtered_bam
    tuple val(sample_id), path("${sample_id}.bwa_se.wrongmap.bam"), emit: ch_se_wrongmap_bam
    tuple val(sample_id), path("${sample_id}.bwa_se.canonical_barcodes.tsv"), emit: ch_se_canonical_barcodes

    script:
    """
    python ${projectDir}/scripts/process_canonical_bam.py --bam_file ${bam} \
                                                          --barcode_file ${barcode} \
                                                          --exon_pos ${exon_pos} \
                                                          --read_type se \
                                                          --soft_clip ${params.filter_softclip_base} \
                                                          --barcode_up ${barcode_up} \
                                                          --barcode_down ${barcode_down} \
                                                          --barcode_check \
                                                          --barcode_temp ${barcode_temp} \
                                                          --output_prefix ${sample_id}.bwa_se \
                                                          --chunk_size 100000 \
                                                          --threads 40
    """
}

process FILTER_PE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(exon_pos), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa_pe.filtered.sorted.bam"), path("${sample_id}.bwa_pe.filtered.sorted.bam.bai"), emit: ch_pe_filtered_bam
    tuple val(sample_id), path("${sample_id}.bwa_pe.wrongmap.r1.fastq.gz"), path("${sample_id}.bwa_pe.wrongmap.r2.fastq.gz"), emit: ch_pe_wrongmap
    tuple val(sample_id), path("${sample_id}.bwa_pe.canonical_barcodes.tsv"), emit: ch_pe_canonical_barcodes
    tuple val(sample_id), path("${sample_id}.bwa_pe.canonical_stats.tsv"), emit: ch_pe_canonical_stats

    script:
    """
    python ${projectDir}/scripts/process_canonical_bam.py --bam_file ${bam} \
                                                          --barcode_file ${barcode} \
                                                          --exon_pos ${exon_pos} \
                                                          --read_type pe \
                                                          --soft_clip ${params.filter_softclip_base} \
                                                          --barcode_up ${barcode_up} \
                                                          --barcode_down ${barcode_down} \
                                                          --barcode_check \
                                                          --barcode_temp ${barcode_temp} \
                                                          --output_prefix ${sample_id}.bwa_pe \
                                                          --chunk_size 100000 \
                                                          --threads 40
    """
}
