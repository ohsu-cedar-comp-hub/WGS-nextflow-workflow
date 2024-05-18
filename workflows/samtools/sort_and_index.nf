#!/usr/bin/env nextflow

include {SortAndIndex} from '../../tools/samtools/sort_and_index.nf'

// Define the workflow for sort and index
workflow { 

    //run sortAndIndex
    SortAndIndex (file(params.bam_unsorted), "test")

}
