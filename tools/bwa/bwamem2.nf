#!/usr/bin/env nextflow

// Define the process for BWA-MEM2 alignment
process bwaMem2Alignment {
    // Set maximum memory
    memory '40 GB'

    // Set output directory for alignment results
    publishDir "${params.outdir}/aligned", mode: 'copy'

    // Define input and output
    input:
    path trim_read1
    path trim_read2
    path idx
    val id 


    output:
    file("${trim_read1.baseName}.bam")

    // BWA-MEM2 alignment command
    script:
    """
    bwa-mem2 mem -K 100000000 -t 6 -Y -M -R "@RG\tID:${id}\tLB:no_library\tPL:illumina\tPU:none\tSM:${trim_read1.baseName}" ${idx} ${trim_read1} ${trim_read2} | samtools view -Sb -@ 4 > ${trim_read1.baseName}.bam
    """
}

