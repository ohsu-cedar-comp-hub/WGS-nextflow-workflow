#!/usr/bin/env nextflow

process RENAME {
    container "${params.container_samtools}"
   
    input:
    path bam

    output:
    path("*")

    script:
    """
    sample_name=\$(samtools view -Sh ${bam} | grep "@RG" | head -1 | cut -f9 | cut -d ':' -f2)
    cp ${bam} \${sample_name}_${bam}
    """
}