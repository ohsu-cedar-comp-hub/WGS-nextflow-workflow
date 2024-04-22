#!/usr/bin/env nextflow

include {GetPileupSummaries} from '../../tools/gatk/get_pileup_summaries.nf'

// Define the workflow for get pileup summaries
workflow {
    
    GetPileupSummaries(file(params.tumor_bam), file(params.normal_bam), file(params.exac), "test")
}
