#!/usr/bin/env nextflow

// filter Mutect2 calls with GATK
process FilterMutectCalls {
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path unfiltered_vcf
    path idx

    output:
    path "${unfiltered_vcf.baseName}_filtered.vcf"

    script:
    """
    gatk FilterMutectCalls -R ${idx} -V ${unfiltered_vcf} -O ${unfiltered_vcf.baseName}_filtered.vcf
    """
}
