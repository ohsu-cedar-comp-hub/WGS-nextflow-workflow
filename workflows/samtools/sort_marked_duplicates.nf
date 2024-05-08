#!/usr/bin/env nextflow

include {SortMarkedDuplicates} from '../../tools/samtools/sort_marked_duplicates.nf'

// Define the workflow for sort and index
workflow { 

    //run sortAndIndex
    SortMarkedDuplicates (file(params.bam_duplicates_unsorted), params.ID)

}