process MATCH_SE_READS {
    label 'process_high_dynamic_memory'
    
    memory {
        def file_size = barcode.size()
        def mem = file_size <= 50_000_000 ? 40 :
                  file_size <= 100_000_000 ? 80 :
                  file_size <= 200_000_000 ? 160 :
                  file_size <= 400_000_000 ? 320 : 640
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(extended_frags), path(ch_hisat2_ref)

    output:
    tuple val(sample_id), path("${sample_id}.match_se.canonical_fail.fastq.gz"), emit: ch_se_canonical_fail
    tuple val(sample_id), path("${sample_id}.match_se.canonical_barcodes.tsv"), emit: ch_se_canonical_barcodes
    tuple val(sample_id), path("${sample_id}.match_se.unknown.tsv"), emit: ch_se_unknown

    script:
    """
    python ${projectDir}/scripts/process_canonical_fastq.py --lib_type ${params.library} \
                                                            --reads ${extended_frags} \
                                                            --read_type se \
                                                            --ref_file ${ch_hisat2_ref} \
                                                            --barcode_file ${barcode} \
                                                            --barcode_up ${barcode_up} \
                                                            --barcode_down ${barcode_down} \
                                                            --barcode_check \
                                                            --barcode_temp ${barcode_temp} \
                                                            --output_prefix ${sample_id}.match_se \
                                                            --chunk_size 100000 \
                                                            --threads 40

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
        py_numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}

process MATCH_PE_READS {
    label 'process_high_dynamic_memory'
    
    memory {
        def file_size = barcode.size()
        def mem = file_size <= 50_000_000 ? 40 :
                  file_size <= 100_000_000 ? 80 :
                  file_size <= 200_000_000 ? 160 :
                  file_size <= 400_000_000 ? 320 : 640
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/canonical_splicing_results/${sample_id}", pattern: "*.canonical_barcodes.tsv", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(not_combined_1), path(not_combined_2), path(ch_hisat2_ref)

    output:
    tuple val(sample_id), path("${sample_id}.match_pe.canonical_fail.r1.fastq.gz"), path("${sample_id}.match_pe.canonical_fail.r2.fastq.gz"), emit: ch_pe_canonical_fail
    tuple val(sample_id), path("${sample_id}.match_pe.canonical_barcodes.tsv"), emit: ch_pe_canonical_barcodes
    tuple val(sample_id), path("${sample_id}.match_se.unknown.tsv"), emit: ch_pe_unknown

    script:
    """
    python ${projectDir}/scripts/process_canonical_fastq.py --lib_type ${params.library} \
                                                            --reads ${extended_frags} \
                                                            --read_type pe \
                                                            --ref_file ${ch_hisat2_ref} \
                                                            --barcode_file ${barcode} \
                                                            --barcode_up ${barcode_up} \
                                                            --barcode_down ${barcode_down} \
                                                            --barcode_check \
                                                            --barcode_temp ${barcode_temp} \
                                                            --output_prefix ${sample_id}.match_pe \
                                                            --chunk_size 100000 \
                                                            --threads 40

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
        py_numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
