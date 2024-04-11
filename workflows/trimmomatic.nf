#!/usr/bin/env nextflow

// define the workflow
workflow {
    // define input paramaters 
    read1=file(params.read1)
    read2=file(params.read2)
    truseq3pefile=file(params.truseq3pefile)

    // run trimmomatic on paired reads
    trimmomaticPE(read1, read2, truseq3pefile)
}
