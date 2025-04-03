/* ---- splicing analysis pipeline ---- */

/* -- load modules -- */
include { check_required } from '../modules/check_software.nf'
include { check_input_files } from '../modules/check_input_files.nf'
include { prepare_files } from '../modules/prepare_files.nf'
include { processing_reads } from '../modules/processing_reads.nf'
include { filtering_reads_se } from '../modules/filtering_reads_se.nf'
include { filtering_reads_pe } from '../modules/filtering_reads_pe.nf'

/* -- define functions -- */
def helpMessage() {
    log.info """
Usage:
    nextflow run nf_splicing/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet        Path of the sample sheet
    
    Optional arguments:
    Basic:
        --outdir              the directory path of output results, default: the current directory
        --library             random, muta, default: muta

    """
}

/* -- initialising parameters -- */
params.help                        = null
params.sample_sheet                = null
params.outdir                      = params.outdir                      ?: "$PWD"
params.library                     = params.library                     ?: "muta"
params.do_pe_reads                 = params.do_pe_reads                 ?: false

params.fastp_cut_mean_quality      = params.fastp_cut_mean_quality      ?: 20
params.flash2_min_overlap          = params.flash2_min_overlap          ?: 10
params.flash2_max_overlap          = params.flash2_max_overlap          ?: 250
params.flash2_min_overlap_outie    = params.flash2_min_overlap_outie    ?: 20
params.flash2_max_mismatch_density = params.flash2_max_mismatch_density ?: 0.25


/* -- check parameters -- */
if (params.help) {
    helpMessage()
    exit 0
}

if (params.sample_sheet) {
    ch_input = Channel.fromPath(file(params.sample_sheet), checkIfExists: true)
                      .splitCsv(header: true, sep: ",")
                      .map { row -> 
                        def sample_id = "${row.sample}_${row.replicate}"
                        tuple(sample_id, row.sample, row.replicate, row.directory, row.read1, row.read2, row.reference, row.barcode) }
} else {
    helpMessage()
    log.info("Error: Please specify the full path of the sample sheet!\n")
    exit 1
}

if (!file(params.outdir).isDirectory()) {
    log.error("Invalid output directory: ${params.outdir}. Please specify a valid directory.")
    exit 1
}

if (params.library != 'random' && params.library != 'muta') {
    log.error("Invalid protocol option: ${params.library}. Valid options: 'random', 'muta'")
    exit 1
}

/* -- check software exist -- */
def required_tools = ['bwa', 'hisat2', 'samtools', 'bamtools', 'flash2', 'fastp']
check_required(required_tools)

/* -- workflow -- */
workflow splicing {
    /* -- check input files exist -- */
    check_input_files(ch_input)
    ch_sample = check_input_files.out.ch_sample_checked

    /* -- prepare the reference files and indexes -- */
    prepare_files(ch_sample)
    ch_bwa_ref = prepare_files.out.ch_bwa_ref
    ch_hisat2_ref = prepare_files.out.ch_hisat2_ref

    /* -- step 1: process reads by fastqc and flash2 -- */
    ch_sample_step1 = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, read1, read2) }
    processing_reads(ch_sample_step1)
    ch_processed_reads = processing_reads.out.ch_merge

    /* -- step 2: align reads to canonical splicing reference by bwa -- */
    ch_sample_step2 = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, barcode) }
                               .join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats ->
                                                                tuple(sample_id,  extended_frags, not_combined_1, not_combined_2) })
                               .join(ch_bwa_ref)

    ch_sample_step2_se = ch_sample_step2.map { sample_id, barcode, extended_frags, not_combined_1, not_combined_2, exon_fasta -> 
                                                tuple(sample_id, barcode, extended_frags, exon_fasta) }
    ch_sample_step2_pe = ch_sample_step2.map { sample_id, barcode, extended_frags, not_combined_1, not_combined_2, exon_fasta -> 
                                                tuple(sample_id, barcode, not_combined_1, not_combined_2, exon_fasta) }

    filtering_reads_se(ch_sample_step2_se)

    if (params.do_pe_reads) {
        filtering_reads_pe(ch_sample_step2_pe)
    }
}
