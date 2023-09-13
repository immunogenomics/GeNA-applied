# Samples SNPs at random across the candidate set, with equal proabability within each MAF decile 

import argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("--snp_info",type=str)
parser.add_argument("--prev_snps",type=str) # to exclude these snps
parser.add_argument("--outfile",type=str)
parser.add_argument("--min_MAF",type=float)
parser.add_argument("--tot_snps",type=int) # must be divisible by 10
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

tot_snps = args.tot_snps
snp_info = args.snp_info
prev_snps = args.prev_snps

snps = pd.read_csv(snp_info, header = None, sep = "\t")
snps.columns = ["ID", "MAF"]

# Ensure the same SNP is not used in both our MAF5 and contMAF sets
prev_snps = pd.read_csv(prev_snps,header = None, sep = "\t").iloc[:,0].values
snps['keep'] = [snps.ID.values[i] not in prev_snps for i in np.arange(snps.shape[0])]
snps = snps.loc[snps.keep,:].reset_index(drop = True)

snps['keep'] = [float(snps.MAF.values[i]) >= args.min_MAF for i in np.arange(snps.shape[0])]
snps = snps.loc[snps.keep,:].reset_index(drop = True)

deciles = np.quantile(snps.MAF.values, np.arange(0.1, 1.1,0.1))
upper_bound = deciles
lower_bound = np.concatenate(([np.min(snps.MAF)],deciles[:9]))

sel_snps = pd.DataFrame({})
for i_decile in np.arange(len(deciles)):
    candidates = snps.loc[snps.MAF.values>=lower_bound[i_decile],:].reset_index(drop = True)
    candidates = candidates.loc[candidates.MAF.values<upper_bound[i_decile],:].reset_index(drop = True)
    i_keep = np.random.choice(candidates.index.values, int(tot_snps/10), replace = False)
    new = pd.DataFrame({"ID": candidates.iloc[i_keep,:]["ID"].values})
    sel_snps = pd.concat((sel_snps, new))
sel_snps.to_csv(args.outfile, header = False, index = False)




