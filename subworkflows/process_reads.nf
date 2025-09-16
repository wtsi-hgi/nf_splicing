workflow process_reads {
    take:
    ch_sample

    main:
    /* -- 1. trim reads -- */
    TRIMMING_READS(ch_sample)
    ch_trim = TRIMMING_READS.out.ch_trim

    /* -- 2. merge reads -- */
    MERGING_READS(ch_trim)
    ch_merge = MERGING_READS.out.ch_merge

    emit:
    ch_merge
}

process TRIMMING_READS {
    label 'process_low'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.r1.fastq.gz"), path("${sample_id}.trimmed.r2.fastq.gz"), path("${sample_id}.trim.txt"), emit: ch_trim

    script:
    """
    fastp --in1               ${read1} \
          --in2               ${read2} \
          --out1              ${sample_id}.trimmed.r1.fastq.gz \
          --out2              ${sample_id}.trimmed.r2.fastq.gz \
          --compression       9 \
          --cut_tail \
          --cut_mean_quality  ${params.fastp_cut_mean_quality} \
          --thread            16 \
          --html              ${sample_id}.trim.html 2>&1 | tee ${sample_id}.trim.txt
    """
}

process MERGING_READS {
    label 'process_low'

    input:
    tuple val(sample_id), path(read1), path(read2), path(trim_stats)

    output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), path("${sample_id}.notCombined_1.fastq.gz"), path("${sample_id}.notCombined_2.fastq.gz"), path("${sample_id}.merge.txt"), path(trim_stats), emit: ch_merge

    script:
    """
    flash2 --min-overlap           ${params.flash2_min_overlap} \
           --max-overlap           ${params.flash2_max_overlap} \
           --min-overlap-outie     ${params.flash2_min_overlap_outie} \
           --max-mismatch-density  ${params.flash2_max_mismatch_density} \
           --output-prefix         ${sample_id} \
           --output-directory      . \
           --threads               32 \
           --allow-outies \
           --compress \
           ${read1} ${read2} 2>&1 | tee ${sample_id}.merge.txt
    """
}
