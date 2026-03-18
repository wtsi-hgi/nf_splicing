process CREATE_SPLICING_COUNTS {
    label 'process_single_dynamic_memory'

    memory {
        def file_size = canonical_barcodes.size()
        def mem = file_size <= 10_000_000 ? 2 :
                  file_size <= 100_000_000 ? 4 :
                  file_size <= 1000_000_000 ? 6 : 8
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/splicing_counts/", pattern: "*.splicing_counts.tsv", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), path(canonical_barcodes), path(classified_junctions)

    output:
    tuple val(sample_id), path("${sample_id}.splicing_counts.tsv"), emit: ch_splicing_counts

    script:
    """
    python ${projectDir}/scripts/create_splicing_counts.py --canonical_file ${canonical_barcodes} \
                                                           --novel_file ${classified_junctions} \
                                                           --barcode_file ${barcode} \
                                                           --output_prefix ${sample_id}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_csv: \$(python -c "import csv; print(csv.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
