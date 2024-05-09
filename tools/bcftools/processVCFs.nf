process processVCFs {

    input:
    path(file)

    output:
    path("${file}_unfiltered_all.vcf.bgz")
    path("${file}_unfiltered_all.vcf.bgz.tbi")

    script:
    """
    files_to_concat=""
    for chr in {1..22} X; do
        files_to_concat+=" ${file}_chr\${chr}_unfiltered.vcf.bgz"
    done

    bcftools view ${file}_unfiltered_all.vcf -Oz -o ${file}_unfiltered_all.vcf.bgz
    bcftools concat -a \$files_to_concat -o ${file}_unfiltered_all.vcf.bgz
    bcftools index -t ${file}_unfiltered_all.vcf.bgz
    """
}
