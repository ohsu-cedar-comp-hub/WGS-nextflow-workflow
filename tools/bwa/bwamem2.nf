#!/usr/bin/env nextflow

// Define the process for BWA-MEM2 alignment
process BWAMEM2 {

    cpus 8
    container "${params.container_bwa}"

    // Set output directory for alignment results
    publishDir "${params.outdir}/aligned/unsorted", mode: 'copy'

    // Define input and output
    input:
    tuple val(sample_id), path(reads)
    path idx
    val id

    output:
    file("${sample_id}.bam")

    // BWA-MEM2 alignment command
    script:
    """
    bwa-mem2 mem -K 100000000 -t ${task.cpus} -Y -M -R "@RG\\tID:${params.id}\\tLB:no_library\\tPL:illumina\\tPU:none\\tSM:${sample_id}" ${params.idx} ${reads[0]} ${reads[1]} | samtools view -Sb -@ ${task.cpus} > ${sample_id}.bam
    """
}

