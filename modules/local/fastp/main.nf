process FASTP {
    label 'process_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.r1.fastq.gz"), path("${sample_id}.trimmed.r2.fastq.gz"), path("${sample_id}.trim.tsv"), emit: ch_trim

    script:
    """
    fastp --in1               ${read1} \
          --in2               ${read2} \
          --out1              ${sample_id}.trimmed.r1.fastq.gz \
          --out2              ${sample_id}.trimmed.r2.fastq.gz \
          --compression       9 \
          --cut_tail \
          --cut_mean_quality  ${params.fastp_cut_mean_quality} \
          --thread            16 \
          --html              ${sample_id}.trim.html 2>&1 | tee ${sample_id}.trim.tsv
    """
}