


include { ANNOTATE } from '../../tools/snpeff/annotate_variants.nf'

workflow {
    // Annotate with snpEff
    ANNOTATE(vcf, sample_id_ch)
}