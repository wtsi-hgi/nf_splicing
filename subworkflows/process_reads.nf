include { FASTP }  from "$projectDir/modules/local/fastp/main"
include { FLASH2 } from "$projectDir/modules/local/flash2/main"

workflow process_reads {
    take:
    ch_sample

    main:
    /* -- 1. trim reads -- */
    FASTP(ch_sample)
    ch_trim = FASTP.out.ch_trim

    /* -- 2. merge reads -- */
    FLASH2(ch_trim)
    ch_merge = FLASH2.out.ch_merge

    emit:
    ch_merge
}
