#!/usr/bin/env nextflow

include { SplitByChromosome } from 
'/Users/goldmael/Desktop/githubWork/WGS-nextflow-workflow/tools/bedtools/split_by_chrom.nf'

// Define the workflow
workflow {
    input_bed=file(params.bed_file)
    SplitByChromosome(input_bed)
}
