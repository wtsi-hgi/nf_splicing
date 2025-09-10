/* ---- splicing analysis pipeline ---- */

/* -- load modules -- */
include { check_required }             from '../modules/check_software.nf'
include { check_input_files }          from '../modules/check_input_files.nf'
include { idxstats_get_values; 
          idxstats_add_values }        from '../modules/format_idxstats.nf'
include { cat_filter_barcodes; 
          cat_map_barcodes;
          cat_beds }                   from '../modules/summary_cat_files.nf'
include { rename_filter_barcodes;
          rename_map_barcodes }        from '../modules/summary_rename_files.nf'
include { hisat2_summary_get_values; 
          hisat2_summary_add_values }  from '../modules/format_hisat2_summary.nf'
include { publish_canonical_barcodes;
          publish_novel_barcodes}      from '../modules/summary_publish_files.nf'

/* -- load subworkflows -- */
include { prepare_files }             from '../subworkflows/prepare_files.nf'
include { process_reads }             from '../subworkflows/process_reads.nf'
include { detect_canonical_se }       from '../subworkflows/detect_canonical_se.nf'
include { detect_canonical_pe }       from '../subworkflows/detect_canonical_pe.nf'
include { detect_novel_se }           from '../subworkflows/detect_novel_se.nf'
include { detect_novel_pe }           from '../subworkflows/detect_novel_pe.nf'
include { summarise_results }         from '../subworkflows/summarise_results.nf'

/* -- define functions -- */
def helpMessage() {
    log.info """
Usage:
    nextflow run nf_splicing/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet                path of the sample sheet
        --library                     random_intron, random_exon, muta_intron, muta_exon, default: random_intron
        --outdir                      the directory path of output results, default: the current directory
    
    Optional arguments:
    Basic:
        --do_pe_reads                 whether to process paired-end reads, default: false
    
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
    """
}

/* -- initialising parameters -- */
params.help                        = null

params.sample_sheet                = null
params.outdir                      = params.outdir                      ?: "$PWD"

params.do_pe_reads                 = params.do_pe_reads                 ?: false
params.library                     = params.library                     ?: "random_intron"

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

params.classify_min_overlap        = params.classify_min_overlap        ?: 2
params.classify_min_cov            = params.classify_min_cov            ?: 2
params.classify_reduce             = params.classify_reduce             ?: 2

/* -- check parameters -- */
if (params.help) {
    helpMessage()
    exit 0
}

if (params.sample_sheet) {
    // reading sample sheet
    def sep = params.sample_sheet.endsWith('.tsv') ? '\t' : ','
    ch_input = Channel.fromPath(file(params.sample_sheet), checkIfExists: true)
                      .splitCsv(header: true, sep: sep)
    
    // check required columns
    def required_cols = ['sample', 'replicate', 'directory', 'read1', 'read2', 'reference', 'barcode', 'barcode_up', 'barcode_down', 'barcode_temp']
    def header_line = new File(params.sample_sheet).readLines().head()
    def header = header_line.split(sep)
    def missing = required_cols.findAll { !(it in header) }

    if (missing) {
        log.error "Sample sheet is missing required columns: ${missing.join(', ')}"
        exit 1
    } else {
        // reformat channel
        ch_input = ch_input.map { row -> 
            def sample_id = "${row.sample}_${row.replicate}"
            tuple(sample_id, row.sample, row.replicate, row.directory, row.read1, row.read2, row.reference, row.barcode, row.barcode_up, row.barcode_down, row.barcode_temp) }
    }
} else {
    helpMessage()
    log.info("Error: Please specify the full path of the sample sheet!\n")
    exit 1
}

def outdir = file(params.outdir)
if (!outdir.exists()) {
    log.info "Output directory does not exist, creating: ${outdir}"
    outdir.mkdirs()
}

if (!file(params.outdir).isDirectory()) {
    log.error("Invalid output directory: ${params.outdir}. Please specify a valid directory.")
    exit 1
}

def valid_library = ['random_intron', 'random_exon', 'muta_intron', 'muta_exon']
if (!(params.library in valid_library)) {
    log.error("Invalid protocol option: ${params.library}. Valid options: ${valid_options.join(', ')}")
    exit 1
}

/* -- check software exist -- */
def required_tools = ['bwa', 'hisat2', 'samtools', 'bamtools', 'flash2', 'fastp']
check_required(required_tools)

/* -- workflow -- */
workflow splicing {
    /* -- check input files exist -- */
    check_input_files(ch_input)
    ch_sample_mapping = check_input_files.out.ch_sample_mapping
    ch_sample_barcodes = check_input_files.out.ch_sample_barcodes

    /* -- prepare the reference files and indexes -- */
    prepare_files(ch_sample_mapping)
    ch_bwa_ref = prepare_files.out.ch_bwa_ref
    ch_hisat2_ref = prepare_files.out.ch_hisat2_ref
    ch_exon_pos = prepare_files.out.ch_exon_pos

    /* -- step 1: process reads by fastp and flash2 -- */
    ch_sample_step1 = ch_sample_mapping.map { sample_id, read1, read2, reference -> tuple(sample_id, read1, read2) }
    process_reads(ch_sample_step1)
    ch_processed_reads = process_reads.out.ch_merge

    /* -- step 2: align reads to canonical splicing reference by bwa -- */
    ch_sample_step2 = ch_sample_barcodes.join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats ->
                                                                        tuple(sample_id, extended_frags, not_combined_1, not_combined_2) })
                                        .join(ch_bwa_ref)
                                        .join(ch_exon_pos)

    ch_sample_step2_se = ch_sample_step2.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos) }
    ch_sample_step2_pe = ch_sample_step2.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos) }

    detect_canonical_se(ch_sample_step2_se)
    ch_bwa_se_fail = detect_canonical_se.out.ch_bwa_se_fail
    ch_bwa_se_filtered_idxstats = detect_canonical_se.out.ch_bwa_se_filtered_idxstats
    ch_bwa_se_barcodes = detect_canonical_se.out.ch_bwa_se_barcodes

    if (params.do_pe_reads) {
        detect_canonical_pe(ch_sample_step2_pe)
        ch_fail_reads_pe = detect_canonical_pe.out.ch_bwa_pe_fail
        ch_bwa_pe_filtered_idxstats = detect_canonical_pe.out.ch_bwa_pe_filtered_idxstats
        ch_bwa_pe_barcodes = detect_canonical_pe.out.ch_bwa_pe_barcodes
    }

    /* -- step 3: align reads to novel splicing reference by hisat2 -- */
    ch_sample_step3_se = ch_sample_step2_se.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos -> 
                                                    tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp) }
                                           .join(ch_bwa_se_fail)
                                           .join(ch_hisat2_ref)
    detect_novel_se(ch_sample_step3_se)
    ch_hisat2_se_summary = detect_novel_se.out.ch_hisat2_se_summary
    ch_hisat2_se_barcodes = detect_novel_se.out.ch_hisat2_se_barcodes
    ch_hisat2_se_fixed = detect_novel_se.out.ch_hisat2_se_fixed
    ch_se_junctions = detect_novel_se.out.ch_se_junctions

    if (params.do_pe_reads) {
        ch_sample_step3_pe = ch_sample_step2_pe.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                        tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp) }
                                               .join(ch_fail_reads_pe)
                                               .join(ch_hisat2_ref)
        detect_novel_pe(ch_sample_step3_pe)
        ch_hisat2_pe_summary = detect_novel_pe.out.ch_hisat2_pe_summary
        ch_hisat2_pe_barcodes = detect_novel_pe.out.ch_hisat2_pe_barcodes
        ch_hisat2_pe_fixed = detect_novel_pe.out.ch_hisat2_pe_fixed
        ch_pe_junctions = detect_novel_pe.out.ch_pe_junctions
    }

    /* -- step 4: summarise results -- */
    
    
    
    // if (params.do_pe_reads) {
    //     idxstats_add_values(ch_bwa_se_filtered_idxstats.join(ch_bwa_pe_filtered_idxstats))
    //     ch_sample_idxstats = idxstats_add_values.out.ch_idxstats

    //     cat_filter_barcodes(ch_bwa_se_barcodes.join(ch_bwa_pe_barcodes))
    //     ch_sample_filter_barcodes = cat_filter_barcodes.out.ch_filter_barcodes

    //     cat_map_barcodes(ch_hisat2_se_barcodes.join(ch_hisat2_pe_barcodes))
    //     ch_sample_map_barcodes = cat_map_barcodes.out.ch_map_barcodes

    //     hisat2_summary_add_values(ch_hisat2_se_summary.join(ch_hisat2_pe_summary))
    //     ch_sample_summary = hisat2_summary_add_values.out.ch_hisat2_summary

    //     cat_beds(ch_se_junctions.join(ch_pe_junctions))
    //     ch_junctions = cat_beds.out.ch_bed
    // } else {
    //     idxstats_get_values(ch_bwa_se_filtered_idxstats)
    //     ch_sample_idxstats = idxstats_get_values.out.ch_idxstats

    //     rename_filter_barcodes(ch_bwa_se_barcodes)
    //     ch_sample_filter_barcodes = rename_filter_barcodes.out.ch_filter_barcodes

    //     rename_map_barcodes(ch_hisat2_se_barcodes)
    //     ch_sample_map_barcodes = rename_map_barcodes.out.ch_map_barcodes

    //     hisat2_summary_get_values(ch_hisat2_se_summary)
    //     ch_sample_summary = hisat2_summary_get_values.out.ch_hisat2_summary

    //     ch_junctions = ch_se_junctions
    // }

    // ch_sample_step4 = ch_input.map { sample_id, sample, replicate, directory, read1, read2, reference, barcode -> tuple(sample_id, sample) }
    //                           .join(ch_sample.map { sample_id, read1, read2, reference, barcode -> tuple(sample_id, barcode) })
    //                           .join(ch_exon_pos)
    //                           .join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats -> 
    //                                                             tuple(sample_id, merge_stats, trim_stats) })   
    //                           .join(ch_sample_idxstats)
    //                           .join(ch_sample_filter_barcodes)
    //                           .join(ch_sample_map_barcodes)
    //                           .join(ch_sample_summary)
    //                           .join(ch_junctions)
                              
    // summarise_results(ch_sample_step4)

    // publish_canonical_barcodes(ch_sample_filter_barcodes)
    // publish_novel_barcodes(ch_sample_map_barcodes)
}
