process CALCULATE_PSI {
    label 'process_single_dynamic_memory'
    
    memory {
        def file_size = splicing_counts[0].size()
        def mem = file_size <= 100_000_000 ? 4 :
                  file_size <= 1_000_000_000 ? 8 :
                  file_size <= 2_000_000_000 ? 16 :
                  file_size <= 4_000_000_000 ? 32 : 64
        "${mem * task.attempt} GB"
    }

    publishDir "${params.outdir}/splicing_reports/${sample}", mode: "copy", overwrite: true

    tag "$sample"

    input:
    tuple val(sample), val(sample_id), val(splicing_counts)

    output:
    tuple val(sample), path("${sample}.canon_only.corrected_psi.tsv"), emit: ch_psi_can_results
    tuple val(sample), path("${sample}.all_events.corrected_psi.tsv"), emit: ch_psi_all_results

    script:
    def list_sample_ids = sample_id.join(',')
    def list_splicing_counts = splicing_counts.join(',')

    """
    ${projectDir}/scripts/calculate_psi_with_error_model.R -r ${projectDir}/scripts \
                                                           -s ${list_sample_ids} \
                                                           -d ${list_splicing_counts} \
                                                           -c can \
                                                           -p ${sample}

    ${projectDir}/scripts/calculate_psi_with_error_model.R -r ${projectDir}/scripts \
                                                           -s ${list_sample_ids} \
                                                           -d ${list_splicing_counts} \
                                                           -c all \
                                                           -p ${sample}

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
