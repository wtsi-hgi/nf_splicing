workflow check_input_files {
    take:
    ch_sample

    main:
    check_files(ch_sample)
    ch_sample_checked = check_files.out.ch_checked

    emit:
    ch_sample_checked
}

process check_files {
    label 'process_single'

    input:
    tuple val(sample_id), val(sample), val(replicate), val(directory), val(read1), val(read2), val(reference), val(barcode)

    output:
    tuple val(sample_id), path("${sample_id}.read_1.fastq.gz"), path("${sample_id}.read_2.fastq.gz"), path("${sample_id}.ref.fasta"), path("${sample_id}.barcodes.txt"), emit: ch_checked

    script:
    def file_read1 = file("${directory}/${read1}")
    def file_read2 = file("${directory}/${read2}")
    def file_reference = file("${directory}/${reference}")
    def file_barcode = file("${directory}/${barcode}")
    
    def valid_read_ext = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
    def valid_ref_ext = [".fa", ".fasta"]

    if (file_read1.exists()) {
        if (!valid_read_ext.any { read1.endsWith(it) }) {
            log.error("Error: File format for ${read1} is incorrect. Expected one of: ${valid_read_ext.join(', ')}")
        }
    } else {
        log.error("Error: ${read1} is not found in ${directory}.")
    }

    if (file_read2.exists()) {
        if (!valid_read_ext.any { read2.endsWith(it) }) {
            log.error("Error: File format for ${read1} is incorrect. Expected one of: ${valid_read_ext.join(', ')}")
        }
    } else {
        log.error("Error: ${read2} is not found in ${directory}.")
    }

    if (file_reference.exists()) {
        if (!valid_ref_ext.any { reference.endsWith(it) }) {
            log.error("Error: File format for ${reference} is incorrect. Expected one of: ${valid_ref_ext.join(', ')}")
        }
    } else {
        log.error("Error: ${reference} is not found in ${directory}.")
    }

    if (file_barcode.exists()) {
        def firstLine = file_barcode.withReader { it.readLine() }
        if (firstLine.contains("\t")) {
            def header = firstLine.split("\t").collect { it.toLowerCase() }
            if (header[0] != "barcode" || header[1] != "varid" || header[2] != "variant" || header[3] != "count") {
                log.error("Error: ${barcode} file format is incorrect. Expected header: barcode\tvarid\tvariant\tcount")
            }
        } else {
            log.error("Error: expect ${barcode} is a tab separted file.")
        }
    } else {
        log.error("Error: ${barcode} is not found in ${directory}.")
    }

    """
    echo "Checking: ${sample_id}"

    if [[ "${file_read1}" == "*.fq" || "${file_read1}" == "*.fastq" ]]
    then
        gzip -c ${file_read1} > ${sample_id}.read_1.fastq.gz
    else
        ln -s ${file_read1} ${sample_id}.read_1.fastq.gz
    fi

    if [[ "${file_read2}" == "*.fq" || "${file_read2}" == "*.fastq" ]]
    then
        gzip -c ${file_read2} > ${sample_id}.read_2.fastq.gz
    else
        ln -s ${file_read2} ${sample_id}.read_2.fastq.gz
    fi

    ln -s ${file_reference} ${sample_id}.ref.fasta
    ln -s ${file_barcode} ${sample_id}.barcodes.txt
    """
}
