workflow processing_reads {
    take:
    ch_sample

    main:
    ch_reads = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, read1, read2) }

    emit:
}

process trimming_reads {
    label 'process_low'

    input:
    tuple val(sample_id), val(read1), val(read2)

    output:
    tuple val(sample_id), path("${sample_id}.r1.fastq.gz"), path("${sample_id}.r2.fastq.gz"), emit: ch_trim

    script:
}

process merging_reads {
    label 'process_low'

    input:
    tuple val(sample_id), val(read1), val(read2)

    output:
    tuple val(sample_id), path("${sample_id}.r1.fastq.gz"), path("${sample_id}.r2.fastq.gz"), emit: ch_merge

    script:
}