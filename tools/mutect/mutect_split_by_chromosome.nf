// Define the process for running MuTect2
process MutectSplitByChromosome {
    // Set maximum memory
    memory '80 G'

    // Set output directory for MuTect2 results
    publishDir "${params.outdir}/svc", mode: 'copy'

    // Define input and output
    input:
    path normal_bam_sorted
    path tumor_bam_sorted
    path chrom
    path mutect_idx
    path pon
    path id
    path germline_resource

    output:
    file("${sample_id}_${chrom}.f1r2.tar.gz")
    file("${sample_id}_${chrom}.vcf")
    file("${sample_id}_${chrom}.vcf.stats"
    // MuTect2 command
    script:
    """
    gatk Mutect2 \
        -R ${params.mutect_idx} \
        -I ${params.normal_bam_sorted} \
        -I ${params.tumor_bam_sorted} \
        -normal "G_${sample_id}" \
        --f1r2-tar-gz "${sample_id}_${chrom}.f1r2.tar.gz" \
        --native-pair-hmm-threads 8 \
        -L "${chrom}.bed" \
        --germline-resource "${germline_resource}" \
        --panel-of-normals ${params.pon} \
        -stats ${sample_id}_${chrom}.vcf.stats \
        -O ${sample_id}_${chrom}.vcf \
    """
}

