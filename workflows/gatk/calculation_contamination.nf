#!/usr/bin/env nextflow

include {CalculateContamination} from '../tools/gatk/calculate_contamination.nf'

// Define the workflow for calculating contamination
workflow {

    bwaMem2Alignment(file(params.tumor_pileup_table), file(params.normal_pileups_table), "test")
}
