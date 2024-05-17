#!/usr/bin/env nextflow

// QC, trimmomatic (optional), alignment with bwa-mem2, and mark duplicates

// import tool modules

include {FastQC} from '../../tools/qc/fastqc/fastqc.nf'
include {TrimmomaticPE} from '../../tools/trimmomatic/trimmomatic.nf'

// define the workflow for fastqc

workflow {
    FastQC(params.read1, params.read2, params.ID)
    TrimmomaticPE(params.read1, params.read2, params.truseq3pefile, params.ID)
}

