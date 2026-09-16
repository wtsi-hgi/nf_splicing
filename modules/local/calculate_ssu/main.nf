process CALCULATE_SSU {
    label 'process_medium'

    publishDir "${params.outdir}/splicing_counts/${sample_id}", mode: "copy", overwrite: true

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

process CORRECT_SSU {
    label 'process_tiny_dynamic_memory'
    
    memory {
        def file_size = ssu_counts[0].size()
        def mem = file_size <= 100_000_000 ? 8 :
                  file_size <= 1_000_000_000 ? 16 :
                  file_size <= 2_000_000_000 ? 32 :
                  file_size <= 4_000_000_000 ? 64 : 128
        "${mem * task.attempt} GB"
    }

    // publishDir "${params.outdir}/splicing_reports/${sample}", mode: "copy", overwrite: true

    tag "$sample"

    input:
    tuple val(sample), val(sample_id), val(ssu_counts)

    output:
    tuple val(sample), path("${sample}.ssu_per_base.details.tsv.gz"), emit: ch_psi_can_results

    script:
    def list_sample_ids = sample_id.join(',')
    def list_ssu_counts = ssu_counts.join(',')    

    """
    ${projectDir}/scripts/calculate_ssu_with_error_model.R -r ${projectDir}/scripts \
                                                           -s ${list_sample_ids} \
                                                           -d ${list_ssu_counts} \
                                                           -p ${sample}

    gzip ${sample}.ssu_per_base.details.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r_base: \$(R --version | grep -i "version" | sed -n '1p' | awk '{print \$3}')
        r_optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))" | tail -n 1 | awk '{print}')
        r_data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" | tail -n 1 | awk '{print}')
        r_glue: \$(Rscript -e "cat(as.character(packageVersion('glue')))" | tail -n 1 | awk '{print}')
        r_tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))" | tail -n 1 | awk '{print}')
        r_vroom: \$(Rscript -e "cat(as.character(packageVersion('vroom')))" | tail -n 1 | awk '{print}')
        r_gtools: \$(Rscript -e "cat(as.character(packageVersion('gtools')))" | tail -n 1 | awk '{print}')
    END_VERSIONS
    """
}
