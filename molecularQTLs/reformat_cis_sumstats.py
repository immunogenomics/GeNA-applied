import pandas as pd
import numpy as np
import argparse
np.random.seed(0)

# Parse Arguments                                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("--infile_hg19_pos",type=str)
parser.add_argument("--infile_sumstats",type=str)
parser.add_argument("--outfile",type=str)

args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

mapped_pos = pd.read_table(args.infile_hg19_pos, header = None)
mapped_pos.columns = ['CHR', 'pos_0', 'POS', 'idx']
sumstats = pd.read_table(args.infile_sumstats)
sumstats = sumstats.iloc[mapped_pos.idx,:] # ignore unmapped loci
sumstats.loc[:,"position"] = mapped_pos.POS.values # update POS to hg19
sumstats = sumstats.loc[:,['chromosome','position','ref','alt', 'maf', 'ma_samples','pvalue', 'beta', 'se', 'gene_id']].drop_duplicates()
sumstats.reset_index(inplace = True, drop = True)

# Rename columns to standards for project
sumstats['CHR'] = sumstats.chromosome.values
sumstats['BP'] = sumstats.position.values
sumstats['P'] = sumstats.pvalue.values
sumstats['BETA'] = sumstats.beta.values
sumstats['SE'] = sumstats.se.values
sumstats['effect_allele'] = sumstats.alt.values
sumstats['other_allele'] = sumstats.ref.values

sumstats = sumstats.loc[:,['CHR','BP','other_allele','effect_allele', 'maf', 'ma_samples','P', 'BETA', 'SE', 'gene_id']]

sumstats.to_csv(args.outfile, index = False)
