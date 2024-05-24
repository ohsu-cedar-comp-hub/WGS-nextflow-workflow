#!/usr/bin/env nextflow

include {CollectAllelicCounts} from '../tools/gatk/collect_allelic_counts.nf'

// Define the workflow for collecting allelic counts
workflow {

    CollectAllelicCounts(file(params.bam_normal_sorted), file(params.idx), file(params.dbsnp), "test")
}

