#!/usr/bin/env nextflow
// process for CollectAllelicCounts


workflow {

    // define input parameters
    bam_normal_sorted = file(params.bam_normal_sorted)
    idx = file(params.idx)
    dbsnp = file(params.dbsnp)

    // run CollectAllelicCounts process
    CollectAllelicCounts(bam_normal_sorted, idx, dbsnp)
}

