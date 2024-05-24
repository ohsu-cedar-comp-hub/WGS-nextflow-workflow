#!/usr/bin/env nextflow

include {MarkDuplicates} from '../../tools/gatk/mark_duplicates.nf'

//Define workflow for marking duplicates
workflow {

    MarkDuplicates(file(params.bam_sorted), params.ID)

}
