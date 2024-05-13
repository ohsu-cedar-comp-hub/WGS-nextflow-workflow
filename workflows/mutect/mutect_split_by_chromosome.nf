#!/usr/bin/env nextflow

include {Mutect2} from '../tools/mutect/mutect_split_by_chromosome.nf'

workflow {
ch_chrom = Channel.from('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X')
ch_stats = Channel.fromPath('*.stats')
ch_f1r2 = Channel.fromPath('*.f1r2.tar.gz')
ch_vcf = Channel.fromPath('*_unfiltered.vcf')

Mutect2(ch_chrom, ch_stats, ch_f1r2, ch_vcf)
}

