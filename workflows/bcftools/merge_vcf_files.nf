#!/usr/bin/env nextflow

include {processVCFs} from '../../tools/bcftools/prepareVCFs.nf
 filesChannel = Channel.fromPath(".vcf")

workflow {
   processVCFs(filesChannel)


}
