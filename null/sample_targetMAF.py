# This script samples SNPs at random from the provided set, including only SNPs with a MAF
# above the specified target MAF, and above the target MAF by no more than 0.005

import argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments                                                                                                 
parser = argparse.ArgumentParser()
parser.add_argument("--infile",type=str) # candidate snps
parser.add_argument("--outfile",type=str)
parser.add_argument("--sample_sz",type=int) # num snps to sample
parser.add_argument("--target_maf",type=float)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

sample_sz = args.sample_sz
infile = args.infile

sel_snps = pd.DataFrame({})
snps = pd.read_csv(infile, header = None, sep = "\t")
snps.columns = ["ID", "MAF"]

MAF = args.target_maf
candidates = snps.loc[snps.MAF.values>MAF,:].reset_index(drop = True)
candidates = candidates.loc[candidates.MAF.values<(MAF+0.005),:].reset_index(drop = True)
i_keep = np.random.choice(candidates.index.values, sample_sz, replace = False)
new = pd.DataFrame({"ID": candidates.iloc[i_keep,:]["ID"].values})
sel_snps = pd.concat((sel_snps, new))
sel_snps.to_csv(args.outfile, header = False, index = False)




