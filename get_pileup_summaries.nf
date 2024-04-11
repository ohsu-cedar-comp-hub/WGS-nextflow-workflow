#!/usr/bin/env nextflow

// Process for getting pileup summaries
process GetPileupSummaries {
    publishDir "${params.outdir}/summaries", mode: 'copy'

    input:
    path tumor_bam
    path exac_vcf

    output:
    file("${tumor_bam.baseName}.pileupsummaries.table")

    script:
    """
    gatk GetPileupSummaries \\
    -I ${params.tumor_bam} \\
    -V ${params.exac_vcf} \\
    -L ${params.exac_vcf} \\
    -O ${tumor_bam.baseName}_pileups_table
    """
}
// Define the workflow
workflow {
    // Define input parameters
    tumor_bam = file(params.tumor_bam)
    exac_vcf = file(params.exac_vcf)

    // Run the GetPileupSummaries process
    GetPileupSummaries(tumor_bam,exac_vcf)
}
