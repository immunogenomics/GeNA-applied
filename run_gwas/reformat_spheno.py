# This script splits the default spheno output for our dataset (with sample IDs formatted as FID_IID)
# to prepare these data as input to plink
import argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("--infile",type=str)
parser.add_argument("--outfile_path",type=str)
args = parser.parse_args()

spheno=pd.read_table(args.infile, index_col=0)
snps = spheno.columns
spheno.insert(0, 'IID',[spheno.index[i].split("_")[1] for i in np.arange(spheno.shape[0])])
spheno.insert(0, 'FID', [spheno.index[i].split("_")[0] for i in np.arange(spheno.shape[0])])
for snp in snps:
    spheno.loc[:,['FID', 'IID', snp]].to_csv(args.outfile_path+snp+".tsv", index=False, sep="\t")
