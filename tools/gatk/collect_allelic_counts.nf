

process CollectAllelicCounts {
    // Set maximum memory
    // memory '40 GB'

    publishDir "${params.outdir}/allelic_counts", mode: 'copy'

    input:
    path bam_normal_sorted
    path idx
    path dbsnp

    output:
    path file("${bam_normal_sorted.baseName}_normal.allelicCounts.csv")

    script:
    """
    gatk CollectAllelicCounts \\
        -L ${dbsnp} \\
        -I ${bam_normal_sorted} \\
        -R ${idx} \\
        -O ${bam_normal_sorted.baseName}_normal.allelicCounts.csv
    """
}
