/* ---- splicing analysis pipeline ---- */

/* -- load modules -- */
include { check_required }    from '../modules/check_software.nf'
include { check_input_files } from '../modules/check_input_files.nf'
include { prepare_files }     from '../modules/prepare_files.nf'
include { process_reads }     from '../modules/process_reads.nf'
include { filter_reads_se }   from '../modules/filter_reads_se.nf'
include { filter_reads_pe }   from '../modules/filter_reads_pe.nf'
include { map_reads_se }      from '../modules/map_reads_se.nf'
include { map_reads_pe }      from '../modules/map_reads_pe.nf'

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
params.do_pe_reads                 = params.do_pe_reads                 ?: false

params.library                     = params.library                     ?: "muta"
params.barcode_template            = params.barcode_template            ?: "NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN"
params.barcode_marker              = params.barcode_marker              ?: "CTACTGATTCGATGCAAGCTTG"

params.fastp_cut_mean_quality      = params.fastp_cut_mean_quality      ?: 20
params.flash2_min_overlap          = params.flash2_min_overlap          ?: 10
params.flash2_max_overlap          = params.flash2_max_overlap          ?: 250
params.flash2_min_overlap_outie    = params.flash2_min_overlap_outie    ?: 20
params.flash2_max_mismatch_density = params.flash2_max_mismatch_density ?: 0.25

params.bwa_gap_open                = params.bwa_gap_open                ?: "10,10"
params.bwa_gap_ext                 = params.bwa_gap_ext                 ?: "5,5"
params.bwa_clip                    = params.bwa_clip                    ?: "1,1"

params.filter_softclip_base        = params.filter_softclip_base        ?: 5

params.hisat2_score_min            = params.hisat2_score_min            ?: "L,0,-0.3"
params.hisat2_mp                   = params.hisat2_mp                   ?: "5,2"
params.hisat2_sp                   = params.hisat2_sp                   ?: "2,1"
params.hisat2_np                   = params.hisat2_np                   ?: 0
params.hisat2_pen_noncansplice     = params.hisat2_pen_noncansplice     ?: 0

params.do_spliced_products         = params.do_spliced_products         ?: false
params.regtools_min_anchor         = params.regtools_min_anchor         ?: 5
params.regtools_min_intron         = params.regtools_min_intron         ?: 20

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
    ch_exon_pos = prepare_files.out.ch_exon_pos

    /* -- step 1: process reads by fastqc and flash2 -- */
    ch_sample_step1 = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, read1, read2) }
    process_reads(ch_sample_step1)
    ch_processed_reads = process_reads.out.ch_merge

    /* -- step 2: align reads to canonical splicing reference by bwa -- */
    ch_sample_step2 = ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, barcode) }
                               .join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats ->
                                                                tuple(sample_id,  extended_frags, not_combined_1, not_combined_2) })
                               .join(ch_bwa_ref)
                               .join(ch_exon_pos)

    ch_sample_step2_se = ch_sample_step2.map { sample_id, barcode, extended_frags, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                tuple(sample_id, barcode, extended_frags, exon_fasta, exon_pos) }
    ch_sample_step2_pe = ch_sample_step2.map { sample_id, barcode, extended_frags, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                tuple(sample_id, barcode, not_combined_1, not_combined_2, exon_fasta, exon_pos) }

    filter_reads_se(ch_sample_step2_se)
    ch_fail_reads_se = filter_reads_se.out.ch_bwa_se_fail

    if (params.do_pe_reads) {
        filter_reads_pe(ch_sample_step2_pe)
        ch_fail_reads_pe = filter_reads_pe.out.ch_bwa_pe_fail
    }

    /* -- step 3: align reads to novel splicing reference by hisat2 -- */
    ch_sample_step3_se = ch_sample_step2_se.map { sample_id, barcode, extended_frags, exon_fasta, exon_pos -> tuple(sample_id, barcode, exon_pos) }
                                           .join(ch_fail_reads_se)
                                           .join(ch_hisat2_ref)
    map_reads_se(ch_sample_step3_se)

    if (params.do_pe_reads) {
        ch_sample_step3_pe = ch_sample_step2_pe.map { sample_id, barcode, not_combined_1, not_combined_2, exon_fasta, exon_pos -> tuple(sample_id, barcode, exon_pos) }
                                               .join(ch_fail_reads_pe)
                                               .join(ch_hisat2_ref)
        map_reads_pe(ch_sample_step3_pe)
    }
}
