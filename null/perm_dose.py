# This script takes an input VCF file and for each SNP generates permutations
# of the genotype dose values across samples

import argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments                                                                                                                                
parser = argparse.ArgumentParser()
parser.add_argument("--infile",type=str)
parser.add_argument("--outfile",type=str)
parser.add_argument("--nperm",type=int)
parser.add_argument("--n_header",type=int)
parser.add_argument("--nampc_file",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

nperm = args.nperm
infile = args.infile
n_header = args.n_header
nampc_file=args.nampc_file

snps = pd.read_csv(infile, skiprows=n_header, sep = "\t")
header = pd.read_csv(infile, sep = "##", header= None, engine='python', nrows=n_header)
n_meta = np.where(snps.columns=='FORMAT')[0][0]+1
meta = snps.iloc[:,:n_meta]
snps = snps.iloc[:,n_meta:]
snps.index=meta.ID.values

# consider only included samples
nampcs=pd.read_csv(nampc_file, index_col=0)
snps = snps.loc[:,nampcs.index]
N = snps.shape[1]

perm_snps = pd.DataFrame({})
for i_snp in np.arange(snps.shape[0]):
    perm_order = np.hstack([np.random.choice(np.arange(N), N, replace=False).reshape(-1,1) for i in np.arange(nperm)])
    G_sub = snps.iloc[i_snp,:].values.reshape(-1,1)
    G_perm = np.hstack([G_sub[perm_order[:,i],:] for i in np.arange(nperm)])
    new = pd.DataFrame(G_perm.T)
    new.index = [snps.index[i_snp].strip()+"_perm"+str(i) for i in np.arange(nperm)]
    perm_snps = pd.concat((perm_snps, new))

perm_snps.columns = snps.columns
perm_meta = pd.DataFrame({'#CHROM': np.repeat(meta['#CHROM'], nperm), 
                          'POS': np.arange(meta.shape[0]*nperm)+100, #dummy, to avoid deletion of matching CHR:POS perm replicates 
                          'ID': perm_snps.index,
                          'REF': np.repeat(meta['REF'], nperm), 
                          'ALT': np.repeat(meta['ALT'], nperm),
                          'QUAL': np.repeat(meta['QUAL'], nperm),
                          'FILTER': np.repeat(meta['FILTER'], nperm),
                          'INFO': np.repeat(meta['INFO'], nperm),
                          'FORMAT': np.repeat(meta['FORMAT'], nperm)})

perm_snps.reset_index(inplace = True, drop = True)
perm_meta.reset_index(inplace = True, drop = True)
export = pd.concat((perm_meta,perm_snps),axis=1)
export.to_csv(args.outfile, header = True, index = False, sep = "\t")
    



