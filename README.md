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
* nextflow
* bwa
* hisat2
* samtools
* bamtools
* R packages
    - Rsamtools


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

# modules
module load HGI/common/nextflow/23.10.0
module load HGI/softpack/users/fs18/nf_longread

#--------------#
# user specify #
#--------------#
# LSF group
export LSB_DEFAULT_USERGROUP=hgi

# Paths
export INPUTSAMPLE=$PWD/inputs/samplesheet.csv
export OUTPUTRES=$PWD/outputs

#-----------#
# pipelines #
#-----------#
nextflow run -resume nf_longread/main.nf --sample_sheet $INPUTSAMPLE \
                                         --protocol DNA \
                                         --platform hifi \
                                         --outdir $OUTPUTRES \
                                         --skip_snvcov
```

<a id="options"></a>

### Usage options
