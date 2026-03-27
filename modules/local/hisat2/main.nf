process HISAT2_ALIGN_SE_READS {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(se_reads), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_se.unique.bam"), emit: ch_se_unique_bam
    tuple val(sample_id), path("${sample_id}.hisat2_se.unmapped.fastq.gz"), emit: ch_se_unmapped
    tuple val(sample_id), path("${sample_id}.hisat2_se.novel_stats.tsv"), emit: ch_se_novel_stats

    script:
    if (params.library in ['muta_exon', 'muta_intron']) {
        """
        python ${projectDir}/scripts/split_fastq_by_ref.py --reads ${se_reads} \
                                                           --read_type se \
                                                           --ref_file ${ref_fasta} \
                                                           --output_prefix ${sample_id} \
                                                           --chunk_size 100000 \
                                                           --threads ${task.cpus}

        for fasta in ${sample_id}.*.exon.fasta; do
            prefix=\${fasta%.exon.fasta}
            fastq=\${prefix}.exon.fastq.gz

            if ! zcat \$fastq 2>/dev/null | head -n 4 | grep -q .; then
                echo "Skipping empty FASTQ: \$fastq"
                continue
            fi

            echo "Processing \$prefix"

            hisat2-build \$fasta \$fasta
            hisat2 -x \$fasta \
                   -U \$fastq \
                   --score-min ${params.hisat2_score_min} \
                   --mp ${params.hisat2_mp} \
                   --sp ${params.hisat2_sp} \
                   --np ${params.hisat2_np} \
                   --pen-noncansplice ${params.hisat2_pen_noncansplice} \
                   --summary-file \${prefix}.hisat2_se.novel_stats.tsv \
                   --new-summary \
                   --threads 32 \
                   -S \${prefix}.hisat2_se.sam

            samtools fastq -@ 32 -f 4 -c 9 -0 \${prefix}.hisat2_se.unmapped.fastq.gz \${prefix}.hisat2_se.sam
            awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==0)||(\$2==16)){print \$0}}' \${prefix}.hisat2_se.sam | grep "NH:i:1\\|^@" | samtools view -b - > \${prefix}.hisat2_se.unique.bam
            rm \${prefix}.hisat2_se.sam
        done

        samtools merge -@ ${task.cpus} ${sample_id}.hisat2_se.unique.bam ${sample_id}.*.hisat2_se.unique.bam
        zcat ${sample_id}.*.hisat2_se.unmapped.fastq.gz > ${sample_id}.hisat2_se.unmapped.fastq.gz

        total=\$(grep "Total reads" *.hisat2_se.novel_stats.tsv | awk -F': ' '{sum += \$2} END{print sum}')
        aligned0=\$(grep "Aligned 0 time" *.hisat2_se.novel_stats.tsv | awk -F': ' '{split(\$2,a," "); sum+=a[1]} END{print sum}')
        aligned1=\$(grep "Aligned 1 time" *.hisat2_se.novel_stats.tsv | awk -F': ' '{split(\$2,a," "); sum+=a[1]} END{print sum}')
        alignedm=\$(grep "Aligned >1 times" *.hisat2_se.novel_stats.tsv | awk -F': ' '{split(\$2,a," "); sum+=a[1]} END{print sum}')

        echo "Total reads: \$total" >> ${sample_id}.hisat2_se.novel_stats.tsv
        echo "Aligned 0 time: \$aligned0" >> ${sample_id}.hisat2_se.novel_stats.tsv
        echo "Aligned 1 time: \$aligned1" >> ${sample_id}.hisat2_se.novel_stats.tsv
        echo "Aligned >1 times: \$alignedm" >> ${sample_id}.hisat2_se.novel_stats.tsv
        echo "Overall alignment rate: \$(( (aligned1 + alignedm)*100 / total ))%" >> ${sample_id}.hisat2_se.novel_stats.tsv

        rm ${sample_id}.*.exon.fasta
        rm ${sample_id}.*.exon.fastq.gz
        rm ${sample_id}.*.hisat2_se.unique.bam
        rm ${sample_id}.*.hisat2_se.unmapped.fastq.gz
        rm ${sample_id}.*.hisat2_se.novel_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
            samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
        END_VERSIONS
        """
    } else {
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
            samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
        END_VERSIONS
        """
    }
}

process HISAT2_ALIGN_PE_READS {
    label 'process_medium'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2), path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.hisat2_pe.unique.bam"), emit: ch_pe_unique_bam
    tuple val(sample_id), path("${sample_id}.hisat2_pe.unmapped.R1.fastq.gz"), path("${sample_id}.hisat2_pe.unmapped.R2.fastq.gz"), emit: ch_pe_unmapped
    tuple val(sample_id), path("${sample_id}.hisat2_pe.novel_stats.tsv"), emit: ch_pe_novel_stats

    script:
    if (params.library in ['muta_exon', 'muta_intron']) {
        """
        python ${projectDir}/scripts/split_fastq_by_ref.py --reads ${read1},${read2} \
                                                           --read_type pe \
                                                           --ref_file ${ref_fasta} \
                                                           --output_prefix ${sample_id} \
                                                           --chunk_size 100000 \
                                                           --threads ${task.cpus}

        for fasta in ${sample_id}.*.exon.fasta; do
            prefix=\${fasta%.exon.fasta}
            fastq1=\${prefix}.exon.R1.fastq.gz
            fastq2=\${prefix}.exon.R2.fastq.gz

            if ! zcat \$fastq1 2>/dev/null | head -n 4 | grep -q .; then
                echo "Skipping empty FASTQ: \$fastq1"
                continue
            fi

            if ! zcat \$fastq2 2>/dev/null | head -n 4 | grep -q .; then
                echo "Skipping empty FASTQ: \$fastq2"
                continue
            fi

            echo "Processing \$prefix"
            hisat2-build \$fasta \$fasta
            hisat2 -x \$fasta \
                   -1 \$fastq1 -2 \$fastq2 --fr \
                   --score-min ${params.hisat2_score_min} \
                   --mp ${params.hisat2_mp} \
                   --sp ${params.hisat2_sp} \
                   --np ${params.hisat2_np} \
                   --pen-noncansplice ${params.hisat2_pen_noncansplice} \
                   --summary-file \${prefix}.hisat2_pe.novel_stats.tsv \
                   --new-summary \
                   --threads 32 \
                   -S \${prefix}.hisat2_pe.sam

            samtools fastq -@ 32 -F 2 -c 9 -1 \${prefix}.hisat2_pe.unmapped.R1.fastq.gz -2 \${prefix}.hisat2_pe.unmapped.R2.fastq.gz -n \${prefix}.hisat2_pe.sam
            awk -F'\\t' -v OFS='\\t' '{if((\$1~/^@/)||(\$2==99)||(\$2==147)||(\$2==83)||(\$2==163)){print \$0}}' \${prefix}.hisat2_pe.sam | grep "NH:i:1\\|^@" | samtools view -b - > \${prefix}.hisat2_pe.unique.bam
            rm \${prefix}.hisat2_pe.sam
        done

        samtools merge -@ ${task.cpus} ${sample_id}.hisat2_pe.unique.bam ${sample_id}.*.hisat2_pe.unique.bam
        zcat ${sample_id}.*.hisat2_pe.unmapped.R1.fastq.gz > ${sample_id}.hisat2_pe.unmapped.R1.fastq.gz
        zcat ${sample_id}.*.hisat2_pe.unmapped.R2.fastq.gz > ${sample_id}.hisat2_pe.unmapped.R2.fastq.gz

        total=\$(grep "Total reads" *.hisat2_pe.novel_stats.tsv | awk -F': ' '{sum += \$2} END{print sum}')
        aligned0=\$(grep "Aligned 0 time" *.hisat2_pe.novel_stats.tsv | awk -F': ' '{split(\$2,a," "); sum+=a[1]} END{print sum}')
        aligned1=\$(grep "Aligned 1 time" *.hisat2_pe.novel_stats.tsv | awk -F': ' '{split(\$2,a," "); sum+=a[1]} END{print sum}')
        alignedm=\$(grep "Aligned >1 times" *.hisat2_pe.novel_stats.tsv | awk -F': ' '{split(\$2,a," "); sum+=a[1]} END{print sum}')

        echo "Total reads: \$total" >> ${sample_id}.hisat2_pe.novel_stats.tsv
        echo "Aligned 0 time: \$aligned0" >> ${sample_id}.hisat2_pe.novel_stats.tsv
        echo "Aligned 1 time: \$aligned1" >> ${sample_id}.hisat2_pe.novel_stats.tsv
        echo "Aligned >1 times: \$alignedm" >> ${sample_id}.hisat2_pe.novel_stats.tsv
        echo "Overall alignment rate: \$(( (aligned1 + alignedm)*100 / total ))%" >> ${sample_id}.hisat2_pe.novel_stats.tsv

        rm ${sample_id}.*.exon.fasta
        rm ${sample_id}.*.exon.R1.fastq.gz
        rm ${sample_id}.*.exon.R2.fastq.gz
        rm ${sample_id}.*.hisat2_pe.unique.bam
        rm ${sample_id}.*.hisat2_pe.unmapped.R1.fastq.gz
        rm ${sample_id}.*.hisat2_pe.unmapped.R2.fastq.gz
        rm ${sample_id}.*.hisat2_pe.novel_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
            samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
        END_VERSIONS
        """
    } else {
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
            samtools: \$( samtools --version | head -n 1 | awk '{print \$2}' )
        END_VERSIONS
        """
    }
}
