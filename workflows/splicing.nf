/* ---- splicing analysis pipeline ---- */

/* -- load modules -- */
include { NOTE_CMD }                  from "$projectDir/modules/local/init_workflow/main"
include { STATS_GET_VALUES; 
          STATS_ADD_VALUES }          from "$projectDir/modules/local/format_stats/main"
include { HISAT2_SUMMARY_GET_VALUES; 
          HISAT2_SUMMARY_ADD_VALUES } from "$projectDir/modules/local/format_hisat2_stats/main"
include { CAT_CANONICAL_BARCODES; 
          CAT_NOVEL_BARCODES;
          CAT_BEDS;
          CAT_BASE_COVS }             from "$projectDir/modules/local/cat_files/main"
include { RENAME_CANONICAL_BARCODES;
          RENAME_NOVEL_BARCODES;
          RENAME_BEDS;
          RENAME_BASE_COVS }          from "$projectDir/modules/local/rename_files/main"

/* -- load subworkflows -- */
include { check_input_files }         from "$projectDir/subworkflows/check_input_files.nf"
include { prepare_files }             from "$projectDir/subworkflows/prepare_files.nf"
include { process_reads }             from "$projectDir/subworkflows/process_reads.nf"
include { detect_canonical_se_align } from "$projectDir/subworkflows/detect_canonical_se_align.nf"
include { detect_canonical_pe_align } from "$projectDir/subworkflows/detect_canonical_pe_align.nf"
include { detect_canonical_se_match } from "$projectDir/subworkflows/detect_canonical_se_match.nf"
include { detect_canonical_pe_match } from "$projectDir/subworkflows/detect_canonical_pe_match.nf"
include { detect_novel_se }           from "$projectDir/subworkflows/detect_novel_se.nf"
include { detect_novel_pe }           from "$projectDir/subworkflows/detect_novel_pe.nf"
include { create_splicing_counts }    from "$projectDir/subworkflows/create_splicing_counts.nf"
include { generate_summary_report }   from "$projectDir/subworkflows/generate_summary_report.nf"

/* -- define functions -- */
def helpMessage() {
    log.info """
Usage:
    nextflow run nf_splicing/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet                path of the sample sheet
        --library                     random_intron, random_exon, random_combi, muta_intron, muta_exon, muta_combi,default: random_intron
        --outdir                      the directory path of output results, default: the current directory
    
    Optional arguments:
    Basic:
        --do_pe_reads                 whether to process paired-end reads, default: false
        --canonical_method            the method of detecting canonical splicing events, align or match, default: match
    
    Fastp:
        --fastp_cut_mean_quality      mean quality for fastp, default: 20
    
    Flash2:
        --flash2_min_overlap          min overlap for flash2, default: 10
        --flash2_max_overlap          max overlap for flash2, default: 250
        --flash2_min_overlap_outie    min overlap outie for flash2, default: 20
        --flash2_max_mismatch_density max mismatch density for flash2, default: 0.25
    
    BWA:
        --bwa_mismatch                mismatch penalty for BWA, default: 4
        --bwa_gap_open                gap open penalty for BWA, default: 10,10
        --bwa_gap_ext                 gap extension penalty for BWA, default: 5,5
        --bwa_clip                    clip penalty for BWA, default: 1,1
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
        --classify_cluster_tol        max tolerance of donor/acceptor positions for clustering, default: 2
        --classify_min_overlap        min anchor to consider partial splicing, default: 2
        --classify_min_cov            min junction coverage to keep, default: 2
    """
}

def check_software_exists(tool) {
    try {
        def process = ["which", tool].execute()
        process.waitFor()
        return process.exitValue() == 0
    } catch (Exception e) {
        return false
    }
}

def check_required(required_tools) {
    def missing_tools = required_tools.findAll { !check_software_exists(it) }

    log.info "====================================="
    log.info "Checking software:"
    required_tools.each { tool ->
        if (check_software_exists(tool)) {
            log.info "    |----> ${tool} is available"
        } else {
            log.info "    |----> ${tool} is not found"
        }
    }

    if (missing_tools) {
        error "Error: the following tools are missing: ${missing_tools.join(', ')}"
    }

    log.info "Done: all required tools are available. Proceeding with the pipeline."
    log.info "====================================="
}

/* -- pipeline info -- */
log.info """
=====================================
${workflow.manifest.name}
Version: ${workflow.manifest.version}
=====================================
"""

/* -- check parameters -- */
if (params.help) {
    helpMessage()
    System.exit(0)
}

if (params.version) {
    println "${workflow.manifest.version}"
    System.exit(0)
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
        error "Error: Sample sheet is missing required columns - ${missing.join(', ')}"
    } else {
        def sheet_file = file(params.sample_sheet)
        log.info("=====================================")
        log.info("Sample sheet content:")
        log.info("-------------------------------------")
        log.info(sheet_file.text)
        log.info("=====================================")

        // reformat channel
        ch_input = ch_input.map { row -> 
            def sample_id = "${row.sample}_${row.replicate}"
            tuple(sample_id, row.sample, row.replicate, row.directory, row.read1, row.read2, row.reference, row.barcode, row.barcode_up, row.barcode_down, row.barcode_temp) }
    }
} else {
    error("Error: Please specify the full path of the sample sheet!\n")
}

def outdir = file(params.outdir)
if (!outdir.exists()) {
    log.info "Output directory does not exist, creating: ${outdir}"
    outdir.mkdirs()
}

if (!file(params.outdir).isDirectory()) {
    error("Invalid output directory: ${params.outdir}. Please specify a valid directory.")
}

def valid_library = ['random_intron', 'random_exon', 'random_combi', 'muta_intron', 'muta_exon', 'muta_combi']
if (!(params.library in valid_library)) {
    error("Invalid library: ${params.library}. Valid options: ${valid_library.join(', ')}")
}

def valid_canonical_method = ['align', 'match']
if (!(params.canonical_method in valid_canonical_method)) {
    error("Invalid canonical method: ${params.canonical_method}. Valid options: ${valid_canonical_method.join(', ')}")
}

/* -- check software exist -- */
def required_tools = ['bwa', 'hisat2', 'samtools', 'bamtools', 'flash2', 'fastp']
check_required(required_tools)

/* -- workflow -- */
workflow splicing {
    /* -- note down the command line -- */
    NOTE_CMD(workflow.commandLine)

    /* -- check input files exist -- */
    check_input_files(ch_input)
    ch_sample_mapping  = check_input_files.out.ch_sample_mapping
    ch_sample_barcodes = check_input_files.out.ch_sample_barcodes

    /* -- prepare the reference files and indexes -- */
    prepare_files(ch_sample_mapping)
    ch_bwa_ref    = prepare_files.out.ch_bwa_ref
    ch_hisat2_ref = prepare_files.out.ch_hisat2_ref
    ch_exon_pos   = prepare_files.out.ch_exon_pos

    /* -- step 1: process reads by fastp and flash2 -- */
    ch_sample_step1 = ch_sample_mapping.map { sample_id, read1, read2, reference -> tuple(sample_id, read1, read2) }
    process_reads(ch_sample_step1)
    ch_processed_reads = process_reads.out.ch_merge

    /* -- step 2: align reads to canonical splicing reference -- */
    if (params.canonical_method == 'match') {
        ch_sample_step2 = ch_sample_barcodes.join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats ->
                                                                        tuple(sample_id, extended_frags, not_combined_1, not_combined_2) })
                                            .join(ch_hisat2_ref)

        ch_sample_step2_se = ch_sample_step2.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, not_combined_1, not_combined_2, ch_hisat2_ref -> 
                                                    tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, ch_hisat2_ref) }
        ch_sample_step2_pe = ch_sample_step2.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, not_combined_1, not_combined_2, ch_hisat2_ref -> 
                                                    tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, ch_hisat2_ref) }

        detect_canonical_se_match(ch_sample_step2_se)
        ch_se_canonical_fail     = detect_canonical_se_match.out.ch_se_canonical_fail
        ch_se_canonical_barcodes = detect_canonical_se_match.out.ch_se_canonical_barcodes
        ch_se_canonical_stats    = detect_canonical_se_match.out.ch_se_canonical_stats

        if (params.do_pe_reads) {
            detect_canonical_pe_match(ch_sample_step2_pe)
            ch_pe_canonical_fail     = detect_canonical_pe_match.out.ch_pe_canonical_fail
            ch_pe_canonical_barcodes = detect_canonical_pe_match.out.ch_pe_canonical_barcodes
            ch_pe_canonical_stats    = detect_canonical_pe_match.out.ch_pe_canonical_stats            
        }
    } else {
        ch_sample_step2 = ch_sample_barcodes.join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats ->
                                                                            tuple(sample_id, extended_frags, not_combined_1, not_combined_2) })
                                            .join(ch_bwa_ref)
                                            .join(ch_exon_pos)

        ch_sample_step2_se = ch_sample_step2.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                    tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, exon_fasta, exon_pos) }
        ch_sample_step2_pe = ch_sample_step2.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp, extended_frags, not_combined_1, not_combined_2, exon_fasta, exon_pos -> 
                                                    tuple(sample_id, barcode, barcode_up, barcode_down, barcode_temp, not_combined_1, not_combined_2, exon_fasta, exon_pos) }

        detect_canonical_se_align(ch_sample_step2_se)
        ch_se_canonical_fail     = detect_canonical_se_align.out.ch_se_canonical_fail
        ch_se_canonical_barcodes = detect_canonical_se_align.out.ch_se_canonical_barcodes
        ch_se_canonical_stats    = detect_canonical_se_align.out.ch_se_canonical_stats

        if (params.do_pe_reads) {
            detect_canonical_pe_align(ch_sample_step2_pe)
            ch_pe_canonical_fail     = detect_canonical_pe_align.out.ch_pe_canonical_fail
            ch_pe_canonical_barcodes = detect_canonical_pe_align.out.ch_pe_canonical_barcodes
            ch_pe_canonical_stats    = detect_canonical_pe_align.out.ch_pe_canonical_stats
        }
    }

    /* -- step 3: align reads to novel splicing reference by hisat2 -- */
    ch_sample_step3_se = ch_sample_barcodes.join(ch_se_canonical_fail)
                                           .join(ch_hisat2_ref)
    detect_novel_se(ch_sample_step3_se)
    ch_se_novel_stats    = detect_novel_se.out.ch_se_novel_stats
    ch_se_novel_barcodes = detect_novel_se.out.ch_se_novel_barcodes
    ch_se_junctions      = detect_novel_se.out.ch_se_junctions
    ch_se_base_cov       = detect_novel_se.out.ch_se_base_cov

    if (params.do_pe_reads) {
        ch_sample_step3_pe = ch_sample_barcodes.join(ch_fail_reads_pe)
                                               .join(ch_hisat2_ref)
        detect_novel_pe(ch_sample_step3_pe)
        ch_pe_novel_stats    = detect_novel_pe.out.ch_pe_novel_stats
        ch_pe_novel_barcodes = detect_novel_pe.out.ch_pe_novel_barcodes
        ch_pe_junctions      = detect_novel_pe.out.ch_pe_junctions
        ch_pe_base_cov       = detect_novel_pe.out.ch_pe_base_cov
    }

    /* -- prepare channels for downstream -- */
    if (params.do_pe_reads) {
        STATS_ADD_VALUES(ch_se_canonical_stats.join(ch_pe_canonical_stats))
        ch_canonical_stats = STATS_ADD_VALUES.out.ch_canonical_stats

        CAT_CANONICAL_BARCODES(ch_se_canonical_barcodes.join(ch_pe_canonical_barcodes))
        ch_canonical_barcodes = CAT_CANONICAL_BARCODES.out.ch_canonical_barcodes

        CAT_NOVEL_BARCODES(ch_se_novel_barcodes.join(ch_pe_novel_barcodes))
        ch_novel_barcodes = CAT_NOVEL_BARCODES.out.ch_novel_barcodes

        HISAT2_SUMMARY_ADD_VALUES(ch_se_novel_stats.join(ch_pe_novel_stats))
        ch_novel_stats = HISAT2_SUMMARY_ADD_VALUES.out.ch_novel_stats

        CAT_BEDS(ch_se_junctions.join(ch_pe_junctions))
        ch_junctions = CAT_BEDS.out.ch_bed

        CAT_BASE_COVS(ch_se_base_cov.join(ch_pe_base_cov))
        ch_base_cov = CAT_BASE_COVS.ch_base_cov
    } else {
        STATS_GET_VALUES(ch_se_canonical_stats)
        ch_canonical_stats = STATS_GET_VALUES.out.ch_canonical_stats

        RENAME_CANONICAL_BARCODES(ch_se_canonical_barcodes)
        ch_canonical_barcodes = RENAME_CANONICAL_BARCODES.out.ch_canonical_barcodes

        RENAME_NOVEL_BARCODES(ch_se_novel_barcodes)
        ch_novel_barcodes = RENAME_NOVEL_BARCODES.out.ch_novel_barcodes

        HISAT2_SUMMARY_GET_VALUES(ch_se_novel_stats)
        ch_novel_stats = HISAT2_SUMMARY_GET_VALUES.out.ch_novel_stats

        RENAME_BEDS(ch_se_junctions)
        ch_junctions = RENAME_BEDS.out.ch_junctions

        RENAME_BASE_COVS(ch_se_base_cov)
        ch_base_cov = RENAME_BASE_COVS.out.ch_base_cov
    }

    /* -- step 4: create count matrices -- */
    ch_sample_step4 = ch_sample_barcodes.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp -> tuple(sample_id, barcode)}
                                        .join(ch_canonical_barcodes)
                                        .join(ch_junctions)
                                        .join(ch_base_cov)
                                        .join(ch_exon_pos)
    create_splicing_counts(ch_sample_step4)
    ch_splicing_counts = create_splicing_counts.out.ch_splicing_counts
    ch_classified_junctions = create_splicing_counts.out.ch_classified_junctions
    ch_splicing_ssu = create_splicing_counts.out.ch_splicing_ssu

    /* -- step 5: summarise results -- */
    ch_sample_step5 = ch_input.map { sample_id, sample, replicate, directory, read1, read2, reference, barcode, barcode_up, barcode_down, barcode_temp -> 
                                        tuple(sample_id, sample) }
                              .join(ch_sample_barcodes.map { sample_id, barcode, barcode_up, barcode_down, barcode_temp ->
                                                                 tuple(sample_id, barcode) })
                              .join(ch_exon_pos)
                              .join(ch_processed_reads.map { sample_id, extended_frags, not_combined_1, not_combined_2, merge_stats, trim_stats -> 
                                                                tuple(sample_id, merge_stats, trim_stats) })   
                              .join(ch_canonical_stats)
                              .join(ch_novel_stats)
                              .join(ch_canonical_barcodes)
                              .join(ch_novel_barcodes)
                              .join(ch_classified_junctions)
                              .join(ch_splicing_counts)

    generate_summary_report(ch_sample_step5)
}
