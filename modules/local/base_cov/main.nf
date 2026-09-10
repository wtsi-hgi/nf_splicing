process BASE_COV {
    label 'process_medium'

    // publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)
    value val(read_type)

    output:
    tuple val(sample_id), path("${sample_id}.base_cov_${read_type}.tsv.gz"), emit: ch_base_cov

    script:
    """
    samtools depth -@ ${task.cpus} -aa ${bam} | pigz -p ${task.cpus} > ${sample_id}.base_cov_${read_type}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}
