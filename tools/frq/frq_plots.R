#!/usr/bin/Rscript

### Example Usage:
# Rscript frq_plots.R -s G_H_507472 -frq vcftools/all_genes.frq -t 2000
# Rscript frq_plots.R -s T_H_50472 -frq vcftools/all_genes.frq -t 2000

# Rscript frq_plots.R -s G_H_507472 -frq vcftools/all_genes.frq -t 1000
# Rscript frq_plots.R -s T_H_50472 -frq vcftools/all_genes.frq -t 1000

# Rscript frq_plots.R -s G_H_507472 -frq vcftools/all_genes.frq -t 4000
# Rscript frq_plots.R -s T_H_50472 -frq vcftools/all_genes.frq -t 4000
###

# TODO update so accepts list of single_sample and loop through "per sample plots" for each sample
  # currently the "all sample plots" are redundant

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-s", "--single_sample", type='character', help="single sample ID")
parser$add_argument("-frq", "--frq_data", type='character', help='allele freq containing all samples')
parser$add_argument('-t', '--threshold',  type='integer', help='threshold of allele freq')
args <- parser$parse_args()

###
# Per sample plots
###
# 1. phred encoded site quality (how much confidence have in variant calls)
var_qual = fread(paste('vcftools/', args$single_sample, '.lqual', sep=''))
a <- ggplot(var_qual, aes(QUAL)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
pdf(paste('sitequal_', args$single_sample,'.pdf',sep=''))
a
dev.off()
# unique(var_qual$QUAL)

# 2. variant mean depth (n reads mapped to each position)
var_depth = fread(paste('vcftools/', args$single_sample, '.ldepth.mean', sep=''))
a <- ggplot(var_depth, aes(MEAN_DEPTH)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
pdf(paste('vardepth_', args$single_sample,'.pdf',sep=''))
a
dev.off()
# tab = summary(var_depth$MEAN_DEPTH)

###
# all sample plots
###
# gene cts for all
df = fread(args$frq_data) %>% as.data.frame()
df = df %>%
  arrange(CHROM, POS)

# freq of all
pdf('gene_freq.pdf')
a = ggplot(df, aes(x=GENE, color=SAMPLE, fill =SAMPLE)) +
  geom_bar() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(title='Gene Frequency (All)', x='Gene', y='Count')
a
dev.off()

# freq of all split by chrom
pdf('gene_freq_byCHROM.pdf')
a+ facet_grid(CHROM ~ .)
dev.off()


# top only (based on combined occurance in all samples)
frq_df = as.data.frame(table(df$GENE) )
colnames(frq_df)= c('GENE', 'Freq')
frq_df = frq_df[order(frq_df$Freq, decreasing = TRUE), ]
row.names(frq_df)= NULL
# filter based on threshold (occurance)
frq_df = frq_df[frq_df$Freq>args$threshold,]
genes = frq_df$GENE %>% as.vector()

# filter orig dataframe for these genes
s1 = subset(df, GENE %in% genes)

# freq of all
pdf(paste('gene_freq_top_',args$threshold,'.pdf',sep=''))
a = ggplot(s1, aes(x=GENE, color=SAMPLE, fill =SAMPLE)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title=paste('Top Gene Frequency (>=',args$threshold, ' Freq)',sep=''), x='Gene', y='Count')
a
dev.off()

# freq of all split by chrom
pdf(paste('gene_freq_byCHROM_top_', args$threshold, '.pdf', sep=''))
a +
  facet_grid(CHROM ~ .) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
