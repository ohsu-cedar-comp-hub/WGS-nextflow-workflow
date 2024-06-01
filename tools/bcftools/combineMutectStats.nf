#!/usr/bin/env nextflow

process MERGESTATS {
    debug true
    
    container "${params.container_gatk}"

    publishDir "${params.outdir}/svc/sort_index", mode: 'copy'
    
    input:
    val stats
    sample_id

    output:
    path "${stats_file.baseName}_all.stats"

    script:
    """
    echo ${stats.join(' ')}
    gatk MergeMutectStats ${stats.join(' ')} -O ${sample_id}_unfiltered.vcf.all.stats
    """
}