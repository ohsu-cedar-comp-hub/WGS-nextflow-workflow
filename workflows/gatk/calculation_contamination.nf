#!/usr/bin/env nextflow

include {CalculateContamination} from '../../tools/gatk/calculate_contamination.nf'

// Define the workflow for calculating contamination
workflow {

    CalculateContamination(file(params.tumor_pileups_table), file(params.normal_pileups_table), "test")
}
