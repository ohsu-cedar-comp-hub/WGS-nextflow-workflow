#!/usr/bin/python

# Example Usage:
# python add_gene.py --sample_ids T_H_50472 G_H_507472 -o vcftools/all_genes.frq

import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--sample_ids', nargs='+', type=str, help='<Required> Set flag, where each sample ID is listed', required=True)
parser.add_argument('-o', '--outfile',required = True, help="output file", type=str)
args = parser.parse_args()

def extract_gene_names(info):
    '''Function to extract gene names from the INFO field. Returns list'''
    gene_names = []
    for item in info.split(';'):
        if item.startswith('ANN='):
            ann_fields = item.split('|')
            gene_name = ann_fields[3]
            gene_names.append(gene_name)
    return gene_names

df_dict = {}

for sample in args.sample_ids :
    vcf_file= 'src/{}.vcf'.format(sample) # original vcf
    file_2 = 'vcftools/{}.frq'.format(sample) # freq table of vcftools

    coord2gene = {}
    # Read VCF file and extract sample names and gene names
    with open(vcf_file, 'r') as f:
        for line in f:

            if line.startswith('#CHROM'):
                sample_names = line.strip().split('\t')[9:]
            elif not line.startswith('#'):
                fields = line.strip().split('\t')
                chrom_col= fields[0]
                pos_col = fields[1]
                info_col = fields[7]
                gene_names_per_sample = extract_gene_names(info_col)
                assert len(gene_names_per_sample)==1, 'Found more than 1 gene for a row'

                # save to look up table
                assert chrom_col+'_'+pos_col not in coord2gene, 'not unique key {} {}'.format(chrom_col, pos_col)
                coord2gene[chrom_col+'_'+pos_col]= gene_names_per_sample

    # add in gene col to output of vcftools
    df2= pd.read_csv(file_2, sep='\t')

    # handle for if ran vcftools with --freq2 then will add a 6th unnamed col
    if df2.index[0].startswith('chr'):
        new_cols = list(df2.columns)+['DETAILS']
        chrom_col = list(df2.index)
        df2.insert(0, 'a', chrom_col)
        df2 = df2.reset_index(drop=True)
        df2.columns= new_cols

    # create list of values for each row
    gene_col = []
    for i in range(0, df2.shape[0]):
        key = df2['CHROM'][i]+'_'+str(df2['POS'][i])
        gene= coord2gene[key][0]
        gene_col.append(gene)
    # add col
    df2['GENE']= gene_col

    # specific to single sample per vcf
    df2['SAMPLE']=[sample]*df2.shape[0]

    df_dict[sample]=df2

# now merge 2 df into one
n = 0
for v in df_dict.values():
    if n==0:
        merged_df = v
        n=1
    else:
        merged_df = pd.concat([merged_df,v],  ignore_index=True)
merged_df.to_csv(args.outfile, sep='\t', index=False)
