/* ---- splicing analysis pipeline ---- */

/* -- load modules -- */
include { check_required } from '../modules/check_software.nf'
include { check_input_files } from '../modules/check_input_files.nf'
include { prepare_files } from '../modules/prepare_files.nf'

/* -- define functions -- */
def helpMessage() {
    log.info """
Usage:
    nextflow run nf_splicing/main.nf --sample_sheet "/path/of/sample/sheet"

    Mandatory arguments:
        --sample_sheet        Path of the sample sheet
        --barcode_file        Path of the barcode association file (relations between barcodes and variants)
    
    Optional arguments:
    Basic:
        --outdir              the directory path of output results, default: the current directory
        --library_type        random, muta, default: muta

    """
}

/* -- check parameters -- */
params.help = null
if (params.help) {
    helpMessage()
    exit 0
}

if (params.sample_sheet) {
    ch_input = Channel.fromPath(file(params.sample_sheet), checkIfExists: true)
                      .splitCsv(header: true, sep: ",")
                      .map { row -> 
                        def sample_id = "${row.sample}_${row.replicate}"
                        tuple(sample_id, row.sample, row.replicate, row.directory, row.read1, row.read2, row.reference) }
} else {
    helpMessage()
    println "Error: Please specify the full path of the sample sheet!\n"
    exit 1
}

params.outdir       = params.outdir       ?: "$PWD"


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
    ch_ref_indexes = prepare_files.out.ch_ref_indexes
    ch_exon_indexes = prepare_files.out.ch_exon_indexes
}
