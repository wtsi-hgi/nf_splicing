process CALCULATE_SSU {
    label 'process_medium'

    publishDir "${params.outdir}/splicing_counts", pattern: "*.splicing_counts.tsv.gz", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(base_cov), path(exon_pos), path(splicing_counts)

    output:
    tuple val(sample_id), path("${sample_id}.splicing_ssu.tsv.gz"), emit: ch_splicing_ssu

    script:
    """
    python ${projectDir}/scripts/calculate_ssu.py --exon_pos ${exon_pos} \
                                                  --splicing_counts ${splicing_counts} \
                                                  --base_cov ${base_cov} \
                                                  --output_dir . \
                                                  --output_prefix ${sample_id} \
                                                  --threads ${task.cpus}
    
    pigz -p ${task.cpus} ${sample_id}.splicing_ssu.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS    
    """
}
