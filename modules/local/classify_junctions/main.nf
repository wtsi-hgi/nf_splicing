process CLASSIFY_NOVEL_JUNCTIONS {
    label 'process_single_dynamic_memory'

    memory {
        def file_size = novel_junctions.size()
        def mem = file_size <= 1_000_000 ? 8 :
                  file_size <= 10_000_000 ? 16 :
                  file_size <= 100_000_000 ? 32 : 64
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.classified_junctions.tsv.gz", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(novel_junctions), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.classified_junctions.tsv.gz"), emit: ch_classified_junctions

    script:
    """
    python ${projectDir}/scripts/classify_novel_junctions.py --lib_type ${params.library} \
                                                             --bed_file ${novel_junctions} \
                                                             --exon_pos ${exon_pos} \
                                                             --cluster_tol ${params.classify_cluster_tol} \
                                                             --junc_cov ${params.classify_min_cov} \
                                                             --min_overlap ${params.classify_min_overlap} \
                                                             --output_prefix ${sample_id}

    gzip ${sample_id}.classified_junctions.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_csv: \$(python -c "import csv; print(csv.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
