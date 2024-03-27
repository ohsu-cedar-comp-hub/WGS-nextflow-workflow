#!/usr/bin/env nextflow

// Define the process for BWA-MEM2 alignment
process bwaMem2Alignment {
    // Set maximum memory
    memory '40 GB'

    // Set output directory for alignment results
    publishDir "${params.outdir}/aligned", mode: 'copy'

    // Define input and output
    input:
    path read1
    path read2
    path idx


    output:
    file("${read1.baseName}.bam")

    // BWA-MEM2 alignment command
    script:
    """
    bwa-mem2 mem -K 1000000 -t 6 -Y -M ${params.idx} ${read1} ${read2} | samtools view -Sb -@ 4 > ${read1.baseName}.bam
    """
}
