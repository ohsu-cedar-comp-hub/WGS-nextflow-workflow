// Merge all f1r2 files into a single tar.gz
process MergeFiles{
    input:
    path vcf_outdir
    path id


    output:
    file("${id}_merged_f1r2.tar.gz")
    file("${id}_merged_stats.txt")
    file("${id}_merged.vcf")

    script:
    """
    //merge f1r2 files
    tar -czf "${id}_merged_f1r2.tar.gz ${vcf_outdir}*f1r2.tar.gz
    //merge stats files
    cat ${vcf_outdir}/svc/*vcf.stats > ${id}_merged_stats.txt"
    //merge vcf files with bcftools
    bcftools concat ${vcf_outdir}*.vcf > ${id}_merged.vcf
    """
}
