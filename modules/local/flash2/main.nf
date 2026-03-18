process FLASH2 {
    label 'process_medium_dynamic_memory'
    
    memory {
        def file_size_1 = read1.size()
        def file_size_2 = read2.size()
        def file_size_total = file_size_1 + file_size_2
        def mem = file_size_total <= 40_000_000_000 ? 20 :
                  file_size_total <= 80_000_000_000 ? 40 :
                  file_size_total <= 160_000_000_000 ? 80 :
                  file_size_total <= 320_000_000_000 ? 160 : 320
        "${mem * task.attempt} GB"
    }

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2), path(trim_stats)

    output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), path("${sample_id}.notCombined_1.fastq.gz"), path("${sample_id}.notCombined_2.fastq.gz"), path("${sample_id}.merge.tsv"), path(trim_stats), emit: ch_merge

    script:
    """
    flash2 --min-overlap           ${params.flash2_min_overlap} \
           --max-overlap           ${params.flash2_max_overlap} \
           --min-overlap-outie     ${params.flash2_min_overlap_outie} \
           --max-mismatch-density  ${params.flash2_max_mismatch_density} \
           --output-prefix         ${sample_id} \
           --output-directory      . \
           --threads               32 \
           --compress \
           ${read1} ${read2} 2>&1 | tee ${sample_id}.merge.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flash2: \$( flash2 --version | head -n 1 | awk '{print \$2}' )
    END_VERSIONS
    """
}