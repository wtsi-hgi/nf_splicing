workflow map_reads_pe {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    hisat2_align_pe_reads(ch_sample)
    ch_hisat2_pe_summary = hisat2_align_pe_reads.out.ch_hisat2_pe_summary
    ch_hisat2_pe_unmapped = hisat2_align_pe_reads.out.ch_hisat2_pe_unmapped
    ch_hisat2_pe_bam = hisat2_align_pe_reads.out.ch_hisat2_pe_bam
    ch_exon_pos = hisat2_align_pe_reads.out.ch_exon_pos

    /* -- 2. fix alignments -- */
    ch_hisat2_pe_bam_join = ch_sample.map { sample_id, barcode, exon_pos, read1, read2, ref_fasta ->
                                            tuple(sample_id, barcode, ref_fasta) }
                                     .join(ch_hisat2_pe_bam)
    fix_pe_reads(ch_hisat2_pe_bam_join)
    ch_hisat2_pe_barcodes = fix_pe_reads.out.ch_hisat2_pe_barcodes
    ch_hisat2_pe_fixed = fix_pe_reads.out.ch_hisat2_pe_fixed

    ch_pe_spliced = params.do_spliced_products 
                            ? fix_pe_reads.out.ch_pe_spliced
                            : Channel.empty()
    
    /* -- 3. extract junctions -- */
    ch_hisat2_pe_fixed_join = ch_hisat2_pe_fixed.join(ch_exon_pos)
    extract_pe_junctions(ch_hisat2_pe_fixed_join)
    ch_pe_junctions = extract_pe_junctions.out.ch_pe_junctions

    emit:
    ch_hisat2_pe_summary
    ch_hisat2_pe_barcodes
    ch_hisat2_pe_fixed
    ch_pe_spliced
    ch_pe_junctions
}

process hisat2_align_pe_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(barcode), path(exon_pos), path(read1), path(read2), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.map_pe.summary.txt"), emit: ch_hisat2_pe_summary
    tuple val(sample_id), path("${sample_id}.map_pe.unmapped.R1.fastq.gz"), path("${sample_id}.map_pe.unmapped.R2.fastq.gz"), emit: ch_hisat2_pe_unmapped
    tuple val(sample_id), path("${sample_id}.map_pe.unique.sorted.bam"), path("${sample_id}.map_pe.unique.sorted.bam.bai"), emit: ch_hisat2_pe_bam
    tuple val(sample_id), path(exon_pos), emit: ch_exon_pos

    script:
    """
    hisat2-build ${ref_fasta} ${ref_fasta}
    hisat2 -x ${ref_fasta} \
           -1 ${read1} -2 ${read2} --fr \
           --score-min ${params.hisat2_score_min} \
           --mp ${params.hisat2_mp} \
           --sp ${params.hisat2_sp} \
           --np ${params.hisat2_np} \
           --pen-noncansplice ${params.hisat2_pen_noncansplice} \
           --summary-file ${sample_id}.map_pe.summary.txt \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.map_pe.sam
    samtools fastq -@ 32 -F 2 -c 9 -1 ${sample_id}.map_pe.unmapped.R1.fastq.gz -2 ${sample_id}.map_pe.unmapped.R2.fastq.gz -n ${sample_id}.map_pe.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==99)||(\$2==147)||(\$2==83)||(\$2==163)){print \$0}}' ${sample_id}.map_pe.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.map_pe.unique.bam
    samtools sort -@ 32 -o ${sample_id}.map_pe.unique.sorted.bam ${sample_id}.map_pe.unique.bam
    samtools index -@ 32 ${sample_id}.map_pe.unique.sorted.bam
    bamtools stats -in ${sample_id}.map_pe.unique.sorted.bam > ${sample_id}.map_pe.unique.txt
    rm ${sample_id}.map_pe.sam ${sample_id}.map_pe.unique.bam
    """
}

process fix_pe_reads {
    label 'process_high'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.bam*", mode: "copy", overwrite: true
    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.spliced_products.txt", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), path(ref_fasta), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.map_pe.barcodes.txt"), emit: ch_hisat2_pe_barcodes
    tuple val(sample_id), path("${sample_id}.map_pe.fixed.sorted.bam"), path("${sample_id}.map_pe.fixed.sorted.bam.bai"), emit: ch_hisat2_pe_fixed
    tuple val(sample_id), path("${sample_id}.map_pe.spliced_products.txt"), emit: ch_pe_spliced

    script:
    def do_spliced_products = params.do_spliced_products ? '-s' : ''

    """
    samtools sort -@ 32 -n -o ${sample_id}.map_pe.unique.sorted_name.bam ${bam}
    ${projectDir}/scripts/fix_bam_pe.R -b ${sample_id}.map_pe.unique.sorted_name.bam -a ${barcode} -r ${ref_fasta} -l ${params.library} ${do_spliced_products}
    mv ${sample_id}.map_pe.unique.sorted_name.barcodes.txt ${sample_id}.map_pe.barcodes.txt
    mv ${sample_id}.map_pe.unique.sorted_name.fixed.sam ${sample_id}.map_pe.fixed.sam

    if [[ -f "${sample_id}.map_pe.unique.sorted_name.spliced_products.txt" ]]; then
        mv ${sample_id}.map_pe.unique.sorted_name.spliced_products.txt ${sample_id}.map_pe.spliced_products.txt
    else
        echo "No data due to missing --do_spliced_products in the pipeline run." > ${sample_id}.map_pe.spliced_products.txt
    fi

    samtools view -@ 64 -b -o ${sample_id}.map_pe.fixed.bam ${sample_id}.map_pe.fixed.sam
    samtools sort -@ 64 -o ${sample_id}.map_pe.fixed.sorted.bam ${sample_id}.map_pe.fixed.bam
    samtools index -@ 64 ${sample_id}.map_pe.fixed.sorted.bam
    rm ${sample_id}.map_pe.fixed.sam ${sample_id}.map_pe.fixed.bam
    """
}

process extract_pe_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_junctions/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.map_pe.junctions.bed"), emit: ch_pe_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.map_pe.junctions.bed ${bam}
    """
}