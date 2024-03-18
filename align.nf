#!/usr/bin/env nextflow

// Define the process for BWA-MEM2 alignment
process bwaMem2Alignment {
    // Set output directory for alignment results
    publishDir "${params.outdir}/aligned", mode: 'copy'

    // Define input and output
    input:
    path reads

    output:
    file("${reads.baseName}.bam")

    // BWA-MEM2 alignment command
    script:
    """
    tar -xvf ${params.ref} && tar -xvf ${params.idx}
    bwa-mem2 mem -t ${task.cpus} ${params.idx}  ${reads} | samtools view -Sb - > ${reads.baseName}.bam
    """
}

// Define the workflow
workflow {
    // Define input parameters
    reads = params.infile

    // Run BWA-MEM2 alignment for each read file
    bwaMem2Alignment(reads)
}
