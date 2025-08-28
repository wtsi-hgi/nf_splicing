workflow prepare_files {
    take:
    ch_sample_mapping

    main:
    ch_ref = ch_sample_mapping.map { sample_id, read1, read2, reference -> tuple(sample_id, reference) }

    create_exon_reference(ch_ref)
    ch_bwa_ref = create_exon_reference.out.ch_bwa_ref
    ch_exon_pos = create_exon_reference.out.ch_exon_pos

    ch_hisat2_ref = ch_ref

    emit:
    ch_bwa_ref
    ch_hisat2_ref
    ch_exon_pos
}

process create_exon_reference {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.exon_ref.fasta"), emit: ch_bwa_ref
    tuple val(sample_id), path("${sample_id}.exon_pos.tsv"), emit: ch_exon_pos

    script:
    """
    python ${params.create_exon_ref_script} -r ${reference} -l ${params.library} -p ${sample_id}
    """
}
