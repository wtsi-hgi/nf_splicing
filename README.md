<div align="center">
<h1 align="center">Splicing Analysis Pipeline</h1>
  <p align="center">A nextflow pipeline of identifying and quantifying splicing events</p>
</div>

## Table of Contents
<details open>
<summary><b>[Show or Hide]</b></summary>

1. [Dependencies](#dependencies)
2. [File Format](#file-format)
    - [Sample sheet](#samplesheet)
3. [Usage](#usage)
    - [Run job](#runjob)
    - [Usage options](#options)
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

    - GenomicRanges = 1.54.1
    - Rsamtools = 2.18.0
    - Biostrings = 2.70.3
</details>

<!-- File Format-->
## File Format

<a id="samplesheet"></a>

### Sample Sheet -- csv
| sample | replicate | directory | read1 | read2 | reference | barcode |
| - | - | - | - | - | - | - |
| s1 | rep1 | /path/of/directory/ | s1_rep1_r1.fastq.gz | s1_rep1_r2.fastq.gz | s1_ref.fa | s1_barcode.txt |
| s1 | rep2 | /path/of/directory/ | s1_rep2_r1.fastq.gz | s1_rep2_r2.fastq.gz | s1_ref.fa | s1_barcode.txt |
| s1 | rep3 | /path/of/directory/ | s1_rep3_r1.fastq.gz | s1_rep3_r2.fastq.gz | s1_ref.fa | s1_barcode.txt |
| s2 | rep1 | /path/of/directory/ | s2_rep1_r1.fastq.gz | s2_rep1_r2.fastq.gz | s2_ref.fa | s2_barcode.txt |
| s2 | rep2 | /path/of/directory/ | s2_rep2_r1.fastq.gz | s2_rep2_r2.fastq.gz | s2_ref.fa | s2_barcode.txt |
| s2 | rep3 | /path/of/directory/ | s2_rep3_r1.fastq.gz | s2_rep3_r2.fastq.gz | s2_ref.fa | s2_barcode.txt |

> [!IMPORTANT]
> 1. The sample sheet must be a csv file and the header must be like above in the example
> 2. all the files should be in the /path/of/directory for each sample


<!-- Usage-->
## Usage

<a id="runjob"></a>

### Run job
submit the bash script below

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

### Usage options
