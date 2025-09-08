workflow map_reads_se {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    hisat2_align_se_reads(ch_sample)
    ch_hisat2_se_summary = hisat2_align_se_reads.out.ch_hisat2_se_summary
    ch_hisat2_se_unmapped = hisat2_align_se_reads.out.ch_hisat2_se_unmapped
    ch_hisat2_se_bam = hisat2_align_se_reads.out.ch_hisat2_se_bam
    ch_exon_pos = hisat2_align_se_reads.out.ch_exon_pos

    /* -- 2. fix alignments -- */
    ch_hisat2_se_bam_join = ch_sample.map { sample_id, barcode, exon_pos, read, ref_fasta ->
                                            tuple(sample_id, barcode, ref_fasta) }
                                     .join(ch_hisat2_se_bam)
    fix_se_reads(ch_hisat2_se_bam_join)
    ch_hisat2_se_barcodes = fix_se_reads.out.ch_hisat2_se_barcodes
    ch_hisat2_se_fixed = fix_se_reads.out.ch_hisat2_se_fixed

    ch_se_spliced = params.do_spliced_products 
                            ? fix_se_reads.out.ch_se_spliced
                            : Channel.empty()
    
    /* -- 3. extract junctions -- */
    ch_hisat2_se_fixed_join = ch_hisat2_se_fixed.join(ch_exon_pos)
    extract_se_junctions(ch_hisat2_se_fixed_join)
    ch_se_junctions = extract_se_junctions.out.ch_se_junctions

    emit:
    ch_hisat2_se_summary
    ch_hisat2_se_barcodes
    ch_hisat2_se_fixed
    ch_se_spliced
    ch_se_junctions
}

process hisat2_align_se_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(barcode), path(exon_pos), path(read), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.summary.txt"), emit: ch_hisat2_se_summary
    tuple val(sample_id), path("${sample_id}.map_se.unmapped.fastq.gz"), emit: ch_hisat2_se_unmapped
    tuple val(sample_id), path("${sample_id}.map_se.unique.sorted.bam"), path("${sample_id}.map_se.unique.sorted.bam.bai"), emit: ch_hisat2_se_bam
    tuple val(sample_id), path(exon_pos), emit: ch_exon_pos

    script:
    """
    hisat2-build ${ref_fasta} ${ref_fasta}
    hisat2 -x ${ref_fasta} \
           -U ${read} \
           --score-min ${params.hisat2_score_min} \
           --mp ${params.hisat2_mp} \
           --sp ${params.hisat2_sp} \
           --np ${params.hisat2_np} \
           --pen-noncansplice ${params.hisat2_pen_noncansplice} \
           --summary-file ${sample_id}.map_se.summary.txt \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.map_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.map_se.unmapped.fastq.gz ${sample_id}.map_se.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==0)||(\$2==16)){print \$0}}' ${sample_id}.map_se.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.map_se.unique.bam
    samtools sort -@ 32 -o ${sample_id}.map_se.unique.sorted.bam ${sample_id}.map_se.unique.bam
    samtools index -@ 32 ${sample_id}.map_se.unique.sorted.bam
    bamtools stats -in ${sample_id}.map_se.unique.sorted.bam > ${sample_id}.map_se.unique.txt
    rm ${sample_id}.map_se.sam ${sample_id}.map_se.unique.bam
    """
}

process fix_se_reads {
    label 'process_high'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.bam*", mode: "copy", overwrite: true
    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.spliced_products.txt", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), path(ref_fasta), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.barcodes.txt"), emit: ch_hisat2_se_barcodes
    tuple val(sample_id), path("${sample_id}.map_se.fixed.sorted.bam"), path("${sample_id}.map_se.fixed.sorted.bam.bai"), emit: ch_hisat2_se_fixed 
    tuple val(sample_id), path("${sample_id}.map_se.spliced_products.txt"), emit: ch_se_spliced
        
    script:
    def do_spliced_products = params.do_spliced_products ? '-s' : ''

    """
    ${projectDir}/scripts/fix_bam_se.R -b ${bam} -a ${barcode} -r ${ref_fasta} -l ${params.library} ${do_spliced_products}
    mv ${sample_id}.map_se.unique.sorted.barcodes.txt ${sample_id}.map_se.barcodes.txt
    mv ${sample_id}.map_se.unique.sorted.fixed.sam ${sample_id}.map_se.fixed.sam

    if [[ -f "${sample_id}.map_se.unique.sorted.spliced_products.txt" ]]; then
        mv ${sample_id}.map_se.unique.sorted.spliced_products.txt ${sample_id}.map_se.spliced_products.txt
    else
        echo "No data due to missing --do_spliced_products in the pipeline run." > ${sample_id}.map_se.spliced_products.txt
    fi

    samtools view -@ 64 -b -o ${sample_id}.map_se.fixed.bam ${sample_id}.map_se.fixed.sam
    samtools sort -@ 64 -o ${sample_id}.map_se.fixed.sorted.bam ${sample_id}.map_se.fixed.bam
    samtools index -@ 64 ${sample_id}.map_se.fixed.sorted.bam
    rm ${sample_id}.map_se.fixed.sam ${sample_id}.map_se.fixed.bam
    """
}

process extract_se_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_junctions/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai), path(exon_pos)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.junctions.bed"), emit: ch_se_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.map_se.junctions.bed ${bam}
    """
}