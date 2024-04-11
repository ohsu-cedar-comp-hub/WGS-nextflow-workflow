#!/usr/bin/env nextflow

workflow {
    // specify input sorted bam
    input_file=file(params.bam_sorted)

    //run MarkDuplications
    MarkDuplicates (input_file)

}
