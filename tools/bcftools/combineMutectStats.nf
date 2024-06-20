#!/usr/bin/env nextflow

process MERGESTATS {
    
    container "${params.container_gatk}"
    
    publishDir "${s3outdir}/svc"
    // publishDir "${params.outdir}/svc/sort_index", mode: 'copy'
    
    input:
    path stats
    val sample_id

    output:
    path "${sample_id}_unfiltered.vcf.all.stats"

    script:
    """
    echo ${stats.join(' ')}
    gatk MergeMutectStats -stats ${stats.join(' -stats ')} -O ${sample_id}_unfiltered.vcf.all.stats
    """
}
