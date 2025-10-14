workflow detect_novel_pe {
    take:
    ch_sample

    main:
    /* -- 1. align reads -- */
    ch_align_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, read1, read2, ref_fasta -> 
                                        tuple(sample_id, read1, read2, ref_fasta) }
    HISAT2_ALIGN_PE_READS(ch_align_input)
    ch_pe_unique_bam = HISAT2_ALIGN_PE_READS.out.ch_pe_unique_bam                               
    ch_pe_novel_stats = HISAT2_ALIGN_PE_READS.out.ch_pe_novel_stats


    /* -- 2. fix alignments -- */
    ch_fix_input = ch_sample.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, read1, read2, ref_fasta -> 
                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, ref_fasta) }
                            .join(ch_pe_unique_bam)
    FIX_PE_READS(ch_fix_input)
    ch_pe_fixed_bam = FIX_PE_READS.out.ch_pe_fixed_bam
    ch_pe_novel_barcodes = FIX_PE_READS.out.ch_pe_novel_barcodes
    
    /* -- 3. extract junctions -- */
    EXTRACT_PE_JUNCTIONS(ch_pe_fixed_bam)
    ch_pe_junctions = EXTRACT_PE_JUNCTIONS.out.ch_pe_junctions

    emit:
    ch_pe_fixed_bam
    ch_pe_novel_barcodes
    ch_pe_junctions
    ch_pe_novel_stats
}

process HISAT2_ALIGN_PE_READS {
    label 'process_medium'

    input:
    tuple val(sample_id), path(read1), path(read2), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.unique.bam"), emit: ch_pe_unique_bam
    tuple val(sample_id), path("${sample_id}.hisat2_pe.unmapped.R1.fastq.gz"), path("${sample_id}.hisat2_pe.unmapped.R2.fastq.gz"), emit: ch_pe_unmapped
    tuple val(sample_id), path("${sample_id}.hisat2_pe.novel_stats.tsv"), emit: ch_pe_novel_stats

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
           --summary-file ${sample_id}.hisat2_pe.novel_stats.tsv \
           --new-summary \
           --threads 32 \
           -S ${sample_id}.hisat2_pe.sam
    samtools fastq -@ 32 -F 2 -c 9 -1 ${sample_id}.hisat2_pe.unmapped.R1.fastq.gz -2 ${sample_id}.hisat2_pe.unmapped.R2.fastq.gz -n ${sample_id}.hisat2_pe.sam
    awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==99)||(\$2==147)||(\$2==83)||(\$2==163)){print \$0}}' ${sample_id}.hisat2_pe.sam | grep "NH:i:1\\|^@" | samtools view -b - > ${sample_id}.hisat2_pe.unique.bam
    rm ${sample_id}.hisat2_pe.sam
    """
}

process FIX_PE_READS {
    label 'process_high_memory'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.bam*", mode: "copy", overwrite: true
    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", pattern: "*.novel_barcodes.tsv", mode: "copy", overwrite: true
 
    input:
    tuple val(sample_id), path(barcode), val(barcode_up), val(barcode_down), val(barcode_temp), path(ref_fasta), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.fixed.sorted.bam"), path("${sample_id}.hisat2_pe.fixed.sorted.bam.bai"), emit: ch_pe_fixed_bam 
    tuple val(sample_id), path("${sample_id}.hisat2_pe.novel_barcodes.tsv"), emit: ch_pe_novel_barcodes

    script:
    def do_spliced_products = params.do_spliced_products ? '--spliced' : ''

    """
    python ${projectDir}/scripts/process_novel_bam.py --lib_type ${params.library} \
                                                      --bam_file ${bam} \
                                                      --ref_file ${ref_fasta} \
                                                      --barcode_file ${barcode} \
                                                      --read_type pe \
                                                      --barcode_up ${barcode_up} \
                                                      --barcode_down ${barcode_down} \
                                                      --barcode_check \
                                                      --barcode_temp ${barcode_temp} \
                                                      ${do_spliced_products} \
                                                      --output_prefix ${sample_id}.hisat2_pe \
                                                      --chunk_size 100000 \
                                                      --threads 40
    samtools sort -@ 64 -o ${sample_id}.hisat2_pe.fixed.sorted.bam ${sample_id}.hisat2_pe.fixed.bam
    samtools index -@ 64 ${sample_id}.hisat2_pe.fixed.sorted.bam
    bamtools stats -in ${sample_id}.hisat2_pe.fixed.sorted.bam > ${sample_id}.hisat2_pe.fixed.tsv
    rm ${sample_id}.hisat2_pe.fixed.bam    
    """
}

process EXTRACT_PE_JUNCTIONS {
    label 'process_single'

    publishDir "${params.outdir}/novel_splicing_results/${sample_id}", mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.junctions.bed"), emit: ch_pe_junctions

    script:
    """
    regtools junctions extract -s RF -a ${params.regtools_min_anchor} -m ${params.regtools_min_intron} -o ${sample_id}.hisat2_pe.junctions.bed ${bam}
    """
}