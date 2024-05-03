#!/usr/bin/env nextflow


// Define the process for running MuTect2
process mutect2 {
    // Set maximum memory
    memory '80 G'

    // Set output directory for MuTect2 results
    publishDir "${params.outdir}/svc", mode: 'copy'

    // Define input and output
    input:
    path normal_bam_sorted
    path tumor_bam_sorted
    path bed_files
    path mutect_idx
    path pon
    path id
    path germline_resource

    output:
    file("${id}_${bed_files}.f1r2.tar.gz")
    file("${id}_${bed_files}.vcf")
    file("${id}_${bed_files}.vcf.stats")

    // MuTect2 command
    script:
    """
    gatk Mutect2 \
        -R ${params.mutect_idx} \
        -I ${params.normal_bam_sorted} \
        -I ${params.tumor_bam_sorted} \
        -normal "G_${id}" \
        --f1r2-tar-gz "${id}_${bed_files}.f1r2.tar.gz" \
        --native-pair-hmm-threads 8 \
        -L "${bed_files}.bed" \
        --germline-resource "${germline_resource}" \
        --panel-of-normals ${pon} \
        -stats ${id}_${bed_files}.vcf.stats
        -O ${id}_${chrom}.vcf \
    """
}
