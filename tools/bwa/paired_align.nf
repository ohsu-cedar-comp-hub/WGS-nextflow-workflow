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
    val id


    output:
    file("${read1.baseName}.bam")

    // BWA-MEM2 alignment command
    script:
    """
    bwa-mem2 mem -K 1000000 -t 6 -Y -M -R "@RG\tID:${id}\tLB:no_library\tPL:illumina\tPU:none\tSM:${read1.baseName}" ${idx} ${read1} ${read2} | samtools view -Sb -@ 4 > ${read1.baseName}.bam
    """
}
