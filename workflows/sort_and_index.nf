#!/usr/bin/env nextflow

workflow { 
    // specify input bam
    input_file=file(params.bam_unsorted)
    //run sortAndIndex
    sortAndIndex (input_file)

}
