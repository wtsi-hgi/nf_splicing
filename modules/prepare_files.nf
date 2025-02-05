workflow prepare_files {
    take:
    ch_sample

    main:
    ch_ref = ch_sample.map { sample_id, sample, replicate, read1, read2, reference -> tuple(sample_id, reference) }

    // build exon indexes
    create_exon_fasta(ch_ref)
    ch_exon = create_exon_fasta.out.ch_exon

    create_bwa_indexes(ch_exon)
    ch_exon_indexes = create_bwa_indexes.out.ch_bwa_indexes

    // build reference indexes
    create_hisat2_indexes(ch_ref)
    ch_ref_indexes = create_hisat2_indexes.out.ch_hisat2_indexes

    emit:
    ch_ref_indexes
    ch_exon_indexes
}

process create_exon_fasta {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.exon.fasta"), emit: ch_exon

    script:
    """
    awk '/^>/ {print \$0; next} {gsub("[acgt]", "", \$0); print \$0}' "${reference}" > "${sample_id}.exon.fasta"
    """
}

process create_hisat2_indexes {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference)

    output:
    tuple val(sample_id), path(reference), path("${sample_id}.ref.*.ht2"), emit: ch_hisat2_indexes

    script:
    """
    hisat2-build "${reference}" "${sample_id}.ref"
    """
}

process create_bwa_indexes {
    label 'process_single'

    input:
    tuple val(sample_id), path(exon_fasta)

    output:
    tuple val(sample_id), path(exon_fasta), path("${exon_fasta}.*"), emit: ch_bwa_indexes

    script:
    """
    bwa index "${exon_fasta}"
    """
}
