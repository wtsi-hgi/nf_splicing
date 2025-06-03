<div align="center">
<h1 align="center">Splicing Analysis Pipeline</h1>
  <p align="center">A nextflow pipeline of identifying and quantifying splicing events</p>
</div>

<p align="center">
  <img src="image/pipeline_structure.gif" alt="Splicing Analysis Pipeline Structure">
</p>

## Table of Contents
<details open>
<summary><b>Catalogue</b></summary>

1. [Dependencies](#dependencies)
2. [File Format](#file-format)
    - [Sample Sheet](#samplesheet)
    - [Reference File](#reference-file)
    - [Barcode File](#barcode-file)
3. [Usage](#usage)
    - [Run](#run)
    - [Options](#options)
4. [Outputs](#outputs)
    - [Structure](#structure)
    - [File Description](#file-description)
5. [Note](#note)
    - [Junction Classification](#junction-classification)
</details>

<!-- Dependencies-->
## Dependencies
<details>
<summary><b>Software</b></summary>

    pigz = 2.7
    BWA = 0.7.17
    HISAT2 = 2.2.1
    Samtools = 1.21
    bamtools = 2.5.2
    FLASH2 = 2.2.00
    fastp = 0.23.4
    RegTools = 1.0.0
    R = 4.3.1
    nextflow = 23.10.0
</details>

<details>
<summary><b>R Packages</b></summary>

    data.table = 1.15.4
    UpSetR = 1.4.0
    gplots = 3.1.3.1
    corrplot = 0.92
    reshape2 = 1.4.4
    optparse = 1.7.4
    psych = 2.4.3
    reactable = 0.4.4
    tidyverse = 2.0.0
    stringr = 1.5.1
    performanceanalytics = 2.0.4
    parallelly = 1.37.1
    dendextend = 1.17.1
    sparkline = 2.0
    ggVennDiagram = 1.5.2
    htmltools = 0.5.8
</details>

<details>
<summary><b>Bioconductor Packages</b></summary>

    GenomicRanges = 1.54.1
    Rsamtools = 2.18.0
    Biostrings = 2.70.3
</details>

<br>

<!-- File Format-->
## File Format

<a id="samplesheet"></a>

### Sample Sheet
| sample | replicate | directory | read1 | read2 | reference | barcode |
| - | - | - | - | - | - | - |
| s1 | rep1 | /path/of/directory/ | s1_rep1_r1.fastq.gz | s1_rep1_r2.fastq.gz | s1_ref.fa | s1_barcode.txt |
| s1 | rep2 | /path/of/directory/ | s1_rep2_r1.fastq.gz | s1_rep2_r2.fastq.gz | s1_ref.fa | s1_barcode.txt |
| s1 | rep3 | /path/of/directory/ | s1_rep3_r1.fastq.gz | s1_rep3_r2.fastq.gz | s1_ref.fa | s1_barcode.txt |
| s2 | rep1 | /path/of/directory/ | s2_rep1_r1.fastq.gz | s2_rep1_r2.fastq.gz | s2_ref.fa | s2_barcode.txt |
| s2 | rep2 | /path/of/directory/ | s2_rep2_r1.fastq.gz | s2_rep2_r2.fastq.gz | s2_ref.fa | s2_barcode.txt |
| s2 | rep3 | /path/of/directory/ | s2_rep3_r1.fastq.gz | s2_rep3_r2.fastq.gz | s2_ref.fa | s2_barcode.txt |

> [!IMPORTANT]
> 1. The sample sheet must be a **csv** file and the header must be like above in the example
> 2. All the files should be in the **directory** for each sample

### Reference File
The reference must be the sequence(s) of the minigene with **exons in uppercase** and **introns in lowercase**.

* **Random Intron Library:** The reference file should contain all variant sequences. For example, if the library includes 100 random intron sequences, the reference fasta file must include all 100 corresponding sequences.
* **Mutagenesis Library:** The reference file only needs the wild-type sequence.

### Barcode File

<br>

<!-- Usage-->
## Usage

<a id="run"></a>

### Run -- Sanger farm
Please submit the bash script below

```bash
#!/bin/bash
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -R "select[mem>1000] rusage[mem=1000]"
#BSUB -M 1000
#BSUB -q normal
#BSUB -J nf_splicing
   
# modules
module load HGI/common/nextflow/23.10.0
module load HGI/softpack/users/fs18/nf_splicing
   
#--------------#
# user specify #
#--------------#
# LSF group
export LSB_DEFAULT_USERGROUP=hgi
   
# Paths
export INPUTSAMPLE=$PWD/inputs/sample_sheet.csv
export OUTPUTRES=$PWD/outputs
  
#-----------#
# pipelines #
#-----------#
nextflow run -resume nf_splicing/main.nf --sample_sheet $INPUTSAMPLE \
                                         --library      random
```

<a id="options"></a>

### Options
#### Mandatory arguments
    --sample_sheet                path of the sample sheet
    --outdir                      the directory path of output results, default: the current directory
    --do_pe_reads                 whether to process paired-end reads, default: false

#### Optional arguments
    Basic:
    --library                     random, muta, default: muta
    --barcode_template            barcode template, default: NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN
    --barcode_marker              barcode marker, default: CTACTGATTCGATGCAAGCTTG

    Fastp:
    --fastp_cut_mean_quality      mean quality for fastp, default: 20
    
    Flash2:
    --flash2_min_overlap          min overlap for flash2, default: 10
    --flash2_max_overlap          max overlap for flash2, default: 250
    --flash2_min_overlap_outie    min overlap outie for flash2, default: 20
    --flash2_max_mismatch_density max mismatch density for flash2, default: 0.25
    
    BWA:
    --bwa_gap_open                gap open penalty for BWA, default: 10,10
    --bwa_gap_ext                 gap extension penalty for BWA, default: 5,5
    --bwa_clip                    clip penalty for BWA, default: 1,1

    Barcode extraction:
    --filter_softclip_base        softclip base for filtering, default: 5
    
    HISAT2:
    --hisat2_score_min            min score for HISAT2, default: L,0,-0.3
    --hisat2_mp                   min/max mismatch penalty for HISAT2, default: 5,2
    --hisat2_sp                   min/max splice penalty for HISAT2, default: 2,1
    --hisat2_np                   non-canonical splicing penalty for HISAT2, default: 0
    --hisat2_pen_noncansplice     non-canonical splicing penalty for HISAT2, default: 0
    
    Spliced products:
    --do_spliced_products         whether to process spliced products, default: false

    Regtools:
    --regtools_min_anchor         min anchor length for regtools, default: 5
    --regtools_min_intron         min intron length for regtools, default: 20

    Junction classification:
    --classify_min_overlap        min overlap for classification, default: 2
    --classify_min_cov            min coverage for classification, default: 2
    --classify_reduce             reduce the number of reads for classification, default: 2


<!-- Outputs-->
## Outputs
### Structure

```bash
📁 output_directory
    ├─── 📁 extracted_barcodes
    │       ├─── 📁 s1_rep1
    │       │       ├─── 📄 canonical_barcodes.txt
    │       │       └─── 📄 novel_barcodes.txt
    │       ├─── 📁 s1_rep2
    │       └─── 📁 s1_rep3
    ├─── 📁 novel_junctions
    │       ├─── 📁 s1_rep1
    │       │       ├─── 📄 junctions.bed
    │       │       ├─── 📄 classified_junctions.txt
    │       │       ├─── 📄 classified_junctions.reduce.txt
    │       │       ├─── 📄 classified_junctions.png
    │       │       └─── 📄 classified_variants.png
    │       ├─── 📁 s1_rep2
    │       └─── 📁 s1_rep3
    ├─── 📁 novel_splicing_results
    │       ├─── 📁 s1_rep1
    │       │       ├─── 📄 spliced_alignment.bam
    │       │       └─── 📄 spliced_products.txt
    │       ├─── 📁 s1_rep2
    │       └─── 📁 s1_rep3
    ├─── 📁 splicing_counts
    │       ├─── 📄 s1_rep1.splicing_matrix.txt
    │       ├─── 📄 s1_rep2.splicing_matrix.txt
    │       └─── 📄 s1_rep3.splicing_matrix.txt
    └─── 📁 splicing_reports
            └─── 📁 s1
                    ├─── 📄 splicing_report.html
                    └─── 📄 *.png
```

### File Description
#### Extracted barcodes
These files summarize all the barcodes in the sequencing library, categorized by canonical and novel splicing events. Compared to the barcode association file, they typically contain more detected barcodes, as the pipeline permits a one-base mismatch during barcode detection.

> [!NOTE]
> 📄 **barcodes of canonical splicing events**
>  
> | sample | replicate | directory | read1 |
> | - | - | - | - |
> | s1 | rep1 | /path/of/directory/ | s1_rep1_r1.fastq.gz |
> | s1 | rep2 | /path/of/directory/ | s1_rep2_r1.fastq.gz |
>
> 📄 **barcodes of novel splicing events**
> 
> | sample | replicate | directory | read1 |
> | - | - | - | - |
> | s1 | rep1 | /path/of/directory/ | s1_rep1_r1.fastq.gz |
> | s1 | rep2 | /path/of/directory/ | s1_rep2_r1.fastq.gz |

<br>

#### Novel junctions

> [!NOTE]
> 📄 **junction bed file**
> 
> This is a bed output file which contains all the junctions detected in the sample.
> 
> 📄 **classified junctions of all the variants**
> 
> This is a tsv output file which contains all the classified junctions for all the variants
> 
> 📄 **classified junctions (reduced)**
>
> This is a tsv output file which contains all the classified junctions regardless of variants, only based on the splicing sites.

<br>

#### Novel splicing results

> [!NOTE]
> 📄 **novel splicing bam**
> 
> This is the bam output file for novel splicing events.
>
> 📄 **novel splicing products**
> 
> This is an output file which contains all the spliced products generating from the bam file.

<br>

#### Splicing counts

This is a raw-read-count file for each splicing event.

<br>

#### Splicing reports

This is a summary html file.

<br>

<!-- Note-->
## Note
### Junction Classification

<p align="center">
  <img src="image/junction_classification.png
  " alt="Junction Classification">
</p>