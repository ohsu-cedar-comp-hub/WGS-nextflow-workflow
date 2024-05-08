#!/usr/bin/env nextflow

include {processVCFs} from '../../tools/bcftools/processVCFs.nf
 filesChannel = Channel.fromPath("/path/to/files.vcf")

workflow {
   processVCFs(filesChannel)


}
