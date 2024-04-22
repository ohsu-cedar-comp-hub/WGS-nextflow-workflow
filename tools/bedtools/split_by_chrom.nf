// define process
process SplitByChromosome {
    publishDir "${params.outdir}/allelic_counts", mode: 'copy'

    input:
    path bed_file

    output:
    path "*.bed"

    script:
    """
    while IFS= read -r line; do
        chr=\$(echo "\$line" | cut -f1)
        if [[ "\$chr" =~ ^(chr[1-9]|chr1[0-9]|chr2[0-2]|chrX)$ ]]; then
            echo "\$line" >> "\${chr}.bed"
        fi
    done < \${bed_file}
    """
}
