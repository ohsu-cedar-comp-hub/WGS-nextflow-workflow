#!/usr/bin/env nextflow

include {MarkDuplicates} from '../tools/gatk/mark_duplicates.nf'

//Define workflow for marking duplicates
workflow {

    MarkDuplicates (file(params.tumor_bam_sorted), file(params.normal_bam_sorted), "test")

}
