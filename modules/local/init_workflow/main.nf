process NOTE_CMD {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val(cmd)

    output:
    path("run_log_*.txt")
    path("sample_sheet_*.tsv")

    script:
    def ts = new Date().format("HHmmss_ddMMMyyyy")
    def module_info = ""
    if (params.sanger_module) {
        module_info = "echo \"Module loaded:\" >> run_log_${ts}.txt"
        module_info += "\nmodule list 2>&1 | head -n -2 >> run_log_${ts}.txt"
    }
    """
    cp ${params.sample_sheet} sample_sheet_${ts}.tsv

    echo "Pipeline: ${params.pipeline_name}" > run_log_${ts}.txt
    echo "Version: ${params.pipeline_version}" >> run_log_${ts}.txt
    echo "Timestamp: ${ts}" >> run_log_${ts}.txt
    echo "" >> run_log_${ts}.txt
    echo "${cmd}" >> run_log_${ts}.txt
    echo "" >> run_log_${ts}.txt
    ${module_info}
    """
}
