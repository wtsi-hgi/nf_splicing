workflow detect_novel_se {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, se_reads, ref_fasta -> 
                                        tuple(sample_id, se_reads, ref_fasta) }
    HISAT2_ALIGN_SE_READS(ch_align_input)
    ch_se_unique_bam = HISAT2_ALIGN_SE_READS.out.ch_se_unique_bam
    ch_se_novel_stats = HISAT2_ALIGN_SE_READS.out.ch_se_novel_stats

    /* -- 2. fix alignments -- */
    ch_fix_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, se_reads, ref_fasta -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, ref_fasta) }
                            .join(ch_se_unique_bam)
    
    FIX_SE_READS(ch_fix_input)
    ch_se_fixed_bam = FIX_SE_READS.out.ch_se_fixed_bam
    ch_se_novel_barcodes = FIX_SE_READS.out.ch_se_novel_barcodes

    SORT_SE_BAM(ch_se_fixed_bam)
    ch_se_sorted_bam = SORT_SE_BAM.out.ch_se_sorted_bam

    /* -- 3. extract junctions -- */
    EXTRACT_SE_JUNCTIONS(ch_se_sorted_bam)
    ch_se_junctions = EXTRACT_SE_JUNCTIONS.out.ch_se_junctions

    emit:
    ch_se_sorted_bam
    ch_se_novel_barcodes
    ch_se_junctions
    ch_se_novel_stats
}

process HISAT2_ALIGN_SE_READS {
    label 'process_medium'

    input:
    tuple val(sample_id), path(se_reads), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.unique.bam"), emit: ch_se_unique_bam
    tuple val(sample_id), path("${sample_id}.hisat2_se.unmapped.fastq.gz"), emit: ch_se_unmapped
    tuple val(sample_id), path("${sample_id}.hisat2_se.novel_stats.tsv"), emit: ch_se_novel_stats

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
           --summary-file ${sample_id}.hisat2_se.novel_stats.tsv \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.hisat2_se.sam
    samtools fastq -@ 32 -f 4 -c 9 -0 ${sample_id}.hisat2_se.unmapped.fastq.gz ${sample_id}.hisat2_se.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==0)||(\$2==16)){print \$0}}' ${sample_id}.hisat2_se.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.hisat2_se.unique.bam
    rm ${sample_id}.hisat2_se.sam
    """
}

process FIX_SE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.novel_barcodes.tsv", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(ref_fasta), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.fixed.bam"), emit: ch_se_fixed_bam 
    tuple val(sample_id), path("${sample_id}.hisat2_se.novel_barcodes.tsv"), emit: ch_se_novel_barcodes

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
                                                      --output_prefix ${sample_id}.hisat2_se \
                                                      --chunk_size 100000 \
                                                      --threads 40
    """
}

process SORT_SE_BAM {
    label 'process_high_memory'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.fixed.sorted.bam"), path("${sample_id}.hisat2_se.fixed.sorted.bam.bai"), emit: ch_se_sorted_bam 

    script:
    """
    samtools sort -@ 40 -o ${sample_id}.hisat2_se.fixed.sorted.bam ${bam}
    samtools index -@ 40 ${sample_id}.hisat2_se.fixed.sorted.bam
    bamtools stats -in ${sample_id}.hisat2_se.fixed.sorted.bam > ${sample_id}.hisat2_se.fixed.tsv
    rm ${bam}
    """
}

process EXTRACT_SE_JUNCTIONS {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.junctions.bed"), emit: ch_se_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.hisat2_se.junctions.bed ${bam}
    """
}