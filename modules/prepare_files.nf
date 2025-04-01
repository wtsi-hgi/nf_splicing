workflow prepare_files {
    take:
    ch_sample

    main:
    ch_ref = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, reference, barcode) }

    /* -- build exon indexes -- */
    create_exon_fasta(ch_ref)
    ch_exon = create_exon_fasta.out.ch_exon

    create_bwa_indexes(ch_exon)
    ch_exon_indexes = create_bwa_indexes.out.ch_bwa_indexes

    /* -- build reference indexes -- */
    if (params.library == 'muta') {
        create_hisat2_indexes_muta(ch_ref)
        ch_ref_indexes = create_hisat2_indexes_muta.out.ch_hisat2_indexes
    } else {
        // random intron needs to re-build the reference
        create_hisat2_indexes_random(ch_ref)
        ch_ref_indexes = create_hisat2_indexes_random.out.ch_hisat2_indexes
    }

    emit:
    ch_ref_indexes
    ch_exon_indexes
}

process create_exon_fasta {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.exon.fasta"), emit: ch_exon

    script:
    """
    awk '/^>/ {print \$0; next} {gsub("[acgt]", "", \$0); print \$0}' "${reference}" > "${sample_id}.exon.fasta"
    """
}

process create_hisat2_indexes_muta {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path(reference), path("${sample_id}.ref.*.ht2"), emit: ch_hisat2_indexes

    script:
    """
    hisat2-build "${reference}" "${sample_id}.ref"
    """
}

process create_hisat2_indexes_random {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.ref.fasta"), path("${sample_id}.ref.*.ht2"), emit: ch_hisat2_indexes

    script:
    """
    awk -F'\\t' -v OFS='\\t' '{if(\$1!="barcode"){print \$2}}' "${barcode}" | sort | uniq > "${sample_id}.random_seq.txt"

    awk -F'\\t' -v OFS='\\t' 'NR==2{exon_count = 0;
                                    intron_count = 0;
                                    while (match(\$0, /[A-Z]+|[a-z]+/)) {
                                        part = substr(\$0, RSTART, RLENGTH);
                                        if (part ~ /^[A-Z]+\$/) {
                                            exon_count++;
                                            exon[exon_count] = part;
                                        } else {
                                            intron_count++;
                                            intron[intron_count] = part;
                                        }
                                        \$0 = substr(\$0, RSTART + RLENGTH);
                                    }
                                    print exon[1]intron[1]exon[2],exon[3]}' "${reference}" > "${sample_id}.constant_seq.txt"

    awk -F'\\t' -v OFS='\\t' 'NR==FNR{part1=\$1;part2=\$2;next}
                                     {print ">var"FNR; print part1\"\"tolower(\$1)\"\"part2}' "${sample_id}.constant_seq.txt" "${sample_id}.random_seq.txt" > ${sample_id}.ref.fasta

    hisat2-build "${sample_id}.ref.fasta" "${sample_id}.ref"
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
