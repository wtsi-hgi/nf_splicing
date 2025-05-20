<div align="center">
<h1 align="center">Splicing Analysis Pipeline</h1>
  <p align="center">A nextflow pipeline of identifying and quantifying splicing events</p>
</div>

## Table of Contents
<details open>
<summary><b>[Show or Hide]</b></summary>

1. [Dependencies](#dependencies)
2. [File Format](#file-format)
    - [Structure of input directories](#structure)
    - [Sample sheet](#samplesheet)
3. [Usage](#usage)
    - [Run job](#runjob)
    - [Usage options](#options)
</details>

<!-- Dependencies-->
## Dependencies
* software
    - pigz = 2.7
    - BWA = 0.7.17
    - HISAT2 = 2.2.1
    - Samtools = 1.21
    - bamtools = 2.5.2
    - FLASH2 = 2.2.00
    - fastp = 0.23.4
    - RegTools = 1.0.0
    - R = 4.3.1
    - nextflow = 23.10.0
* R packages
    - data.table = 1.15.4
    - UpSetR = 1.4.0
    - gplots = 3.1.3.1
    - corrplot = 0.92
    - reshape2 = 1.4.4
    - optparse = 1.7.4
    - psych = 2.4.3
    - reactable = 0.4.4
    - tidyverse = 2.0.0
    - stringr = 1.5.1
    - performanceanalytics = 2.0.4
    - parallelly = 1.37.1
    - dendextend = 1.17.1
    - sparkline = 2.0
    - ggVennDiagram = 1.5.2
    - htmltools = 0.5.8
* Bioconductor packages
    - GenomicRanges = 1.54.1
    - Rsamtools = 2.18.0
    - Biostrings = 2.70.3
    - 

<!-- File Format-->
## File Format

<a id="structure"></a>

### Structure of input directories
![example](./image/inputs.png)

<a id="samplesheet"></a>

### Sample Sheet -- csv
| sample | replicate | directory | read1 | read2 | reference |
| - | - | - | - | - | - |
| s1 | rep1 | /path/of/directory/ | s1rep1_r1.fastq.gz | s1rep1_r2.fastq.gz | ref.fa |
| s1 | rep2 | /path/of/directory/ | s1rep2_r1.fastq.gz | s1rep2_r2.fastq.gz | ref.fa |
| s1 | rep3 | /path/of/directory/ | s1rep3_r1.fastq.gz | s1rep3_r2.fastq.gz | ref.fa |

> [!Note]  
> 1. The sample sheet must be a csv file and the header must be like below in the example


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
