#!/usr/bin/env nextflow


workflow {
    // specify input sorted bam
    input_file=file(params.svc_vcf)

    //run MarkDuplications
    Annotate_Variants (input_file)

}
