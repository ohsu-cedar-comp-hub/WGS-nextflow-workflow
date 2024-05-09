#!/usr/bin/env nextflow

process combineStatsFiles {
    input:
    val chromosome
    path stats_file

    output:
    path "${stats_file.baseName}_all.stats"

    script:
    """
    stats_options=""
    for stat_file in "${stats_file.baseName}"_chr"${chromosome}"_unfiltered.vcf.stats; do
        stats_options+="-stats \$stat_file "
    done

    echo "Running command: gatk MergeMutectStats \$stats_options -O ${stats_file.baseName}_all.stats"

    gatk MergeMutectStats \$stats_options -O ${stats_file.baseName}_all.stats
    """
}
