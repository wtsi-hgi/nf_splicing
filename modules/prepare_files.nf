workflow prepare_files {
    take:
    ch_sample

    main:
    ch_ref = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, reference, barcode) }

    /* -- create bwa reference of exon fasta -- */
    create_bwa_reference(ch_ref)
    ch_bwa_ref = create_bwa_reference.out.ch_bwa_ref

    /* -- create hisat2 reference of wildtype or random intron -- */
    if (params.library == 'muta') {
        ch_hisat2_ref = ch_ref.map { sample_id, reference, barcode -> tuple(sample_id, reference) }
    } else {
        // random intron needs to re-build the reference
        create_hisat2_reference(ch_ref)
        ch_hisat2_ref = create_hisat2_reference.out.ch_hisat2_ref
    }

    emit:
    ch_bwa_ref
    ch_hisat2_ref
}

process create_bwa_reference {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.exon.fasta"), emit: ch_bwa_ref

    script:
    """
    awk '/^>/ {print \$0; next} {gsub("[acgt]", "", \$0); print \$0}' "${reference}" > "${sample_id}.exon.fasta"
    """
}

process create_hisat2_reference {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.ref.fasta"), emit: ch_hisat2_ref

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
    """
}
