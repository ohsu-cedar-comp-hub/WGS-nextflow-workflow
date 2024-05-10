#!/usr/bin/env nextflow

workflow {
    include { prepareVCFs } from '../../tools/bcftools/prepareVCFs.nf'
    Channel
        .fromPath("{params.outdir}*.vcf") // unknown if this will work as desired
        .set { unfiltered_vcfs }
        
    unfiltered_vcfs
        .map { vcf -> tuple(vcf) }
        .into { prepare_input }

    prepareVCFs(prepare_input)

    Channel
        .of((1..22).collect { it.toString() } + 'X')
        .set { chromosomes }
