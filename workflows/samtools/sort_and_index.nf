#!/usr/bin/env nextflow

include {sortAndIndex} from '../tools/samtools/sort_and_index.nf'

// Define the workflow for sort and index
workflow { 

    //run sortAndIndex
    sortAndIndex (file(params.bam_unsorted), "test")

}
