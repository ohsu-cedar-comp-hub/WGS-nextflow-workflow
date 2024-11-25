#!/usr/bin/env nextflow

process GETINTERVALS {
    container "${params.container_samtools}"

    input:
    path bam
    path bai

    output:
    stdout

    script:
    """
    samtools view -SH ${bam} | grep "@SQ" | cut -d ":" -f2 | cut -f1 | sed -n '/^\\(chr1\\|chr2\\|chr3\\|chr4\\|chr5\\|chr6\\|chr7\\|chr8\\|chr9\\|chr10\\|chr11\\|chr12\\|chr13\\|chr14\\|chr15\\|chr16\\|chr17\\|chr18\\|chr19\\|chr20\\|chr21\\|chr22\\|chr23\\|chr24\\|chrX\\|chrY\\|chrUn\\|chrN\\)/p' | sed '/decoy\$/d'
    """
}

