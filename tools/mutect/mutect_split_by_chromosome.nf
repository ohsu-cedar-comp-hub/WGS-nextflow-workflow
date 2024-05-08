#!/usr/bin/env nextflow


// Define the process for running MuTect2
process mutect2 {
    // Set maximum memory
    memory '80 G'

    // Set output directory for MuTect2 results
    publishDir "${params.outdir}/svc", mode: 'copy'

    // Define input and output
   input:
    path id
    path mutect_idx
    path chrom
    path tumor_bam_sorted
    path normal_bam_sorted

    output:
    file "${id}_${chrom}_unfiltered.vcf"
    file "${id}_${chrom}_f1r2.tar.gz"
    file "${id}_${chrom}_unfiltered.vcf.stats"

    script:
    """
    gatk mutect2 \\
        -R ${mutect_idx} \\
        -I ${tumor_bam_sorted} \\
        -I ${normal_bam_sorted} \\
        --panel-of-normals ${pon} \\
        -normal ${normal_bam.baseName} \\
        -L ${chrom} \\
        -O ${id}_${chrom}_unfiltered.vcf \\
        --f1r2-tar-gz ${id}_${chrom}_f1r2.tar.gz \\
        -stats ${id}_${chrom}_unfiltered.vcf.stats
    """
}
