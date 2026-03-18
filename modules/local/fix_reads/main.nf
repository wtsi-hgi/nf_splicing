process FIX_SE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.novel_barcodes.tsv", mode: "copy", overwrite: true

    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(ref_fasta), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.fixed.bam"), emit: ch_se_fixed_bam 
    tuple val(sample_id), path("${sample_id}.hisat2_se.novel_barcodes.tsv"), emit: ch_se_novel_barcodes

    script:
    def do_spliced_products = params.do_spliced_products ? '--spliced' : ''

    """
    python ${projectDir}/scripts/process_novel_bam.py --lib_type ${params.library} \
                                                      --bam_file ${bam} \
                                                      --ref_file ${ref_fasta} \
                                                      --barcode_file ${barcode} \
                                                      --read_type se \
                                                      --barcode_up ${barcode_up} \
                                                      --barcode_down ${barcode_down} \
                                                      --barcode_check \
                                                      --barcode_temp ${barcode_temp} \
                                                      ${do_spliced_products} \
                                                      --output_prefix ${sample_id}.hisat2_se \
                                                      --chunk_size 100000 \
                                                      --threads 40

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_csv: \$(python -c "import csv; print(csv.__version__)")
        py_pysam: \$(python -c "import pysam; print(pysam.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
        py_pandas: \$(python -c "import pandas; print(pandas.__version__)")
        py_Bio: \$(python -c "import Bio; print(Bio.__version__)")
        py_numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}

process FIX_PE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.novel_barcodes.tsv", mode: "copy", overwrite: true
 
    tag "$sample_id"

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(ref_fasta), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.fixed.bam"), emit: ch_pe_fixed_bam 
    tuple val(sample_id), path("${sample_id}.hisat2_pe.novel_barcodes.tsv"), emit: ch_pe_novel_barcodes

    script:
    def do_spliced_products = params.do_spliced_products ? '--spliced' : ''

    """
    python ${projectDir}/scripts/process_novel_bam.py --lib_type ${params.library} \
                                                      --bam_file ${bam} \
                                                      --ref_file ${ref_fasta} \
                                                      --barcode_file ${barcode} \
                                                      --read_type pe \
                                                      --barcode_up ${barcode_up} \
                                                      --barcode_down ${barcode_down} \
                                                      --barcode_check \
                                                      --barcode_temp ${barcode_temp} \
                                                      ${do_spliced_products} \
                                                      --output_prefix ${sample_id}.hisat2_pe \
                                                      --chunk_size 100000 \
                                                      --threads 40
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        py_argparse: \$(python -c "import argparse; print(argparse.__version__)")
        py_csv: \$(python -c "import csv; print(csv.__version__)")
        py_pysam: \$(python -c "import pysam; print(pysam.__version__)")
        py_polars: \$(python -c "import polars; print(polars.__version__)")
        py_pandas: \$(python -c "import pandas; print(pandas.__version__)")
        py_Bio: \$(python -c "import Bio; print(Bio.__version__)")
        py_numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
