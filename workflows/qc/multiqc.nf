#!/usr/bin/env nextflow

include {MultiQC} from '../../tools/qc/multiqc/multiqc.nf'

// Define the workflow for multiqc
workflow {

    MultiQC(file(params.fastqc_read1), file(params.fastqc_read2), file(params.trim_fastqc_read1), file(params.trim_fastqc_read2), params.ID)
}
