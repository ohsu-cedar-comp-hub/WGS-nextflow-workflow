#!/usr/bin/bash

# Get allele freq
bash allele_frq.sh src/G_H_507472.vcf vcftools/G_H_507472
bash allele_frq.sh src/T_H_50472.vcf vcftools/T_H_50472

# Add gene symbol info
python add_gene.py --sample_ids T_H_50472 G_H_507472 -o vcftools/all_genes.frq

# Generate plots
Rscript frq_plots.R -s G_H_507472 -frq vcftools/all_genes.frq -t 2000
Rscript frq_plots.R -s T_H_50472 -frq vcftools/all_genes.frq -t 2000

Rscript frq_plots.R -s G_H_507472 -frq vcftools/all_genes.frq -t 1000
Rscript frq_plots.R -s T_H_50472 -frq vcftools/all_genes.frq -t 1000

Rscript frq_plots.R -s G_H_507472 -frq vcftools/all_genes.frq -t 4000
Rscript frq_plots.R -s T_H_50472 -frq vcftools/all_genes.frq -t 4000
