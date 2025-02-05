#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { splicing } from './workflows/splicing.nf'

workflow {
    splicing()
}
