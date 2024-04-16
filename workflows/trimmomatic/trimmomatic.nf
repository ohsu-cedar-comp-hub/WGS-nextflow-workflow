#!/usr/bin/env nextflow

include {trimmomaticPE} from '../tools/trimmomatic/trimmomatic.nf'

// define the workflow for trimmomaticPE
workflow {
   
    trimmomaticPE (file(params.read1), file(params.read2), file(params.truseq3pefile), "test")
}
