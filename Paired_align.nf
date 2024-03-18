#!/usr/bin/env nextflow

// Define the process for BWA-MEM2 alignment
process bwaMem2Alignment {
    // Set output directory for alignment results
    publishDir "${params.outdir}/aligned", mode: 'copy'

    // Define input and output
    input:
    path Read1
    path Read2

    output:
    file("${Read1.baseName}.bam")

    // BWA-MEM2 alignment command
    script:
    """
    tar -xvf ${params.idx}
    bwa-mem2 mem -t ${task.cpus} ${params.idx}  ${Read1} ${Read2} | samtools view -Sb - > ${Read1.baseName}.bam
    """
}

// Define the workflow
workflow {
    // Define input parameters
    Read1=file(params.Read1)
    Read2=file(params.Read2)

    // Run BWA-MEM2 alignment for each read file
    bwaMem2Alignment(Read1, Read2)
}
