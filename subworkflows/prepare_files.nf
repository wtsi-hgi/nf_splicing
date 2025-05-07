workflow prepare_files {
    take:
    ch_sample

    main:
    ch_ref = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, reference, barcode) }

    /* -- 1. create bwa reference of exon fasta -- */
    create_bwa_reference(ch_ref)
    ch_bwa_ref = create_bwa_reference.out.ch_bwa_ref

    /* -- 2. user has to create hisat2 reference by himself -- */
    ch_hisat2_ref = ch_ref.map { sample_id, reference, barcode -> tuple(sample_id, reference) }

    /* -- 3. create exon positions -- */
    create_exon_positions(ch_ref)
    ch_exon_pos = create_exon_positions.out.ch_exon_pos

    emit:
    ch_bwa_ref
    ch_hisat2_ref
    ch_exon_pos
}

process create_bwa_reference {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.exon.fasta"), emit: ch_bwa_ref

    script:
    """
    awk -F'\\t' -v OFS='\\t' 'NR==2{
        split(\$0, gene_part, "NNNN");
        exon_count = 0;
        exon_seq = "";
        exon[1] = ""; exon[2] = ""; exon[3] = "";
        for (i = 1; i <= length(gene_part[1]); i++) {
            char = substr(gene_part[1], i, 1);
            if (char ~ /[A-Z]/) {
                exon_seq = exon_seq char;
            } else if (exon_seq != "") {
                exon_count++;
                exon[exon_count] = exon_seq;
                exon_seq = "";
            }
        }
        if (exon_seq != "") {
            exon_count++;
            exon[exon_count] = exon_seq;
        }
        print ">E1_E2_E3";
        print exon[1] exon[2] exon[3];
        print ">E1_E3";
        print exon[1] exon[3];
    }' "${reference}" > "${sample_id}.exon.fasta"
    """
}

process create_exon_positions {
    label 'process_single'

    input:
    tuple val(sample_id), path(reference), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.exon_pos.txt"), emit: ch_exon_pos

    script:
    """
    awk -F'\\t' -v OFS='\\t' 'NR==2{
        split(\$0, gene_part, "NNNN");
        exon_start = 0;
        exon_end = 0;
        is_exon = 0;
        j = 1;
        for (i = 1; i <= length(gene_part[1]); i++) {
            char = substr(gene_part[1], i, 1);
            if (char ~ /[A-Z]/) {
                if (!is_exon) {
                    exon_start = i;
                    is_exon = 1;
                }
                exon_end = i;
            } else if (is_exon) {
                print "E"j, exon_start, exon_end;
                is_exon = 0;
                j++;
            }
        }
        if (is_exon) {
            print "E"j, exon_start, exon_end;
            j++;
        }
    }' "${reference}" > "${sample_id}.exon_pos.txt"
    """
}