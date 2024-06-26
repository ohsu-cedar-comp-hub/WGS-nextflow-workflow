#!/usr/bin/env nextflow

include {TrimmomaticPE} from '../../tools/trimmomatic/trimmomatic.nf'

// define the workflow for trimmomaticPE
workflow {
   
    TrimmomaticPE (file(params.read1), file(params.read2), file(params.truseq3pefile), "test")
}
