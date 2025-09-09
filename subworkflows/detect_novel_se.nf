workflow detect_novel_se {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, se_reads, ref_fasta -> 
                                        tuple(sample_id, se_reads, ref_fasta) }
    hisat2_align_se_reads(ch_align_input)
    ch_hisat2_se_summary = hisat2_align_se_reads.out.ch_hisat2_se_summary
    ch_hisat2_se_unmapped = hisat2_align_se_reads.out.ch_hisat2_se_unmapped
    ch_hisat2_se_bam = hisat2_align_se_reads.out.ch_hisat2_se_bam

    /* -- 2. fix alignments -- */
    ch_fix_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, se_reads, ref_fasta -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, ref_fasta) }
                            .join(ch_hisat2_se_bam)
    fix_se_reads(ch_fix_input)
    ch_hisat2_se_barcodes = fix_se_reads.out.ch_hisat2_se_barcodes
    ch_hisat2_se_fixed = fix_se_reads.out.ch_hisat2_se_fixed
    
    /* -- 3. extract junctions -- */
    extract_se_junctions(ch_hisat2_se_fixed)
    ch_se_junctions = extract_se_junctions.out.ch_se_junctions

    emit:
    ch_hisat2_se_summary
    ch_hisat2_se_barcodes
    ch_hisat2_se_fixed
    ch_se_junctions
}

process hisat2_align_se_reads {
    label 'process_medium'

    input:
    tuple val(sample_id), path(se_reads), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.summary.txt"), emit: ch_hisat2_se_summary
    tuple val(sample_id), path("${sample_id}.map_se.unmapped.fastq.gz"), emit: ch_hisat2_se_unmapped
    tuple val(sample_id), path("${sample_id}.map_se.unique.bam"), emit: ch_hisat2_se_bam

    script:
    """
    hisat2-build ${ref_fasta} ${ref_fasta}
    hisat2 -x ${ref_fasta} \
           -U ${se_reads} \
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
    rm ${sample_id}.map_se.sam
    """
}

process fix_se_reads {
    label 'process_high_memory'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.barcodes.txt", mode: "copy", overwrite: true
    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.bam*", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(ref_fasta), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.barcodes.txt"), emit: ch_hisat2_se_barcodes
    tuple val(sample_id), path("${sample_id}.map_se.fixed.sorted.bam"), path("${sample_id}.map_se.fixed.sorted.bam.bai"), emit: ch_hisat2_se_fixed 
  
    script:
    def do_spliced_products = params.do_spliced_products ? '--spliced' : ''

    """
    python ${projectDir}/scripts/process_novel_bam.py --lib_type ${params.library} \
                                                      --bam_file ${bam} \
                                                      --ref_file ${ref_fasta} \
                                                      --barcode_file ${barcode} \
                                                      --read_type se \
                                                      --barcode_up ${barcode_up} \
                                                      --barcode_down ${barcode_down} \
                                                      --barcode_check \
                                                      --barcode_temp ${barcode_temp} \
                                                      ${do_spliced_products} \
                                                      --output_prefix ${sample_id}.map_se \
                                                      --chunk_size 100000 \
                                                      --threads 40
    samtools sort -@ 64 -o ${sample_id}.map_se.fixed.sorted.bam ${sample_id}.map_se.fixed.bam
    samtools index -@ 64 ${sample_id}.map_se.fixed.sorted.bam
    bamtools stats -in ${sample_id}.map_se.fixed.sorted.bam > ${sample_id}.map_se.fixed.txt
    rm ${sample_id}.map_se.fixed.bam    
    """
}

process extract_se_junctions {
    label 'process_single'

    publishDir "${params.outdir}/novel_junctions/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.map_se.junctions.bed"), emit: ch_se_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.map_se.junctions.bed ${bam}
    """
}