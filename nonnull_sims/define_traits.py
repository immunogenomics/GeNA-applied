# This script defines trait values per individual for several real phenotypes in the given 
# single-cell data object. 
# 
# Global gene expression program traits: top [n_hPCs] gene expression PCs
# 
# Cluster abundance traits: fractional abundance of all cluster-based cell types
#       Only clusters uncorrelated with batch (Pearson’s r2 < 0.25 to any batch)
#       Only clusters with representation of at least 50 cells from at least 100 samples
#
# Cell type specific gene expression programs: top [n_celltype_PCs] gene expression PCs among cells in the cluster
#       Only for clusters that passed the QC metrics for cluster abundance traits

import argparse
import numpy as np
import pandas as pd
import cna
import scipy.stats as st
import anndata as ad
import scanpy as sc
np.random.seed(0)

# Parse Arguments                                                                                                                                                
parser = argparse.ArgumentParser()
parser.add_argument("--majortype",type=str)
parser.add_argument("--sc_object",type=str)
parser.add_argument("--outfile",type=str)
parser.add_argument("--n_hPCs",type=int,default=10)
parser.add_argument("--n_celltype_PCs",type=int,default=3)
parser.add_argument("--clustQC_minsamples",type=int,default=100)
parser.add_argument("--clustQC_mincells",type=str,default=50)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

n_hPCs = args.n_hPCs
n_celltype_PCs = args.n_celltype_PCs

# import real single-cell data object
d = cna.read(args.sc_object)

if args.majortype=="allcells":
    celltype_col = 'majortype'
else:
    celltype_col = 'celltype'

# cluster QC
celltypes = d.obs[celltype_col].value_counts().index
for celltype in celltypes:
    d.obs["count_"+celltype] = 1*(d.obs[celltype_col].values==celltype)
trait_candidates = ["count_"+celltypes[i] for i in np.arange(len(celltypes))]
d.obs_to_sample(trait_candidates, aggregate = np.sum)
keep = np.sum(d.samplem[trait_candidates]>args.clustQC_mincells, axis=0).values>args.clustQC_minsamples
trait_candidates = np.array(trait_candidates)[keep]
celltypes = [trait_candidates[i].split("count_")[1] for i in np.arange(len(trait_candidates))]
print("Retained celltype clusters:")
print(celltypes)

# cell type fractional abundance                                                                                                                              
for celltype in celltypes:
    d.obs["frac_"+celltype] = 1*(d.obs[celltype_col].values==celltype)
traits = ["frac_"+celltypes[i] for i in np.arange(len(celltypes))]
d.obs_to_sample(traits, aggregate = np.mean)

# remove any cluster with r-sq>0.25 to any batch
onehot_batches = pd.get_dummies(d.samplem['batch'])
batch_check = d.samplem[traits].join(onehot_batches)
batch_corrs = batch_check.corr()[traits].iloc[len(traits):,:]
keep = np.sum(batch_corrs**2>0.25, axis=0)==0
traits = np.array(traits)[keep.values].tolist()
celltypes = np.array(celltypes)[keep.values].tolist()

# gene expression program
for i_hPC in np.arange(n_hPCs):
    d.obs['hPC'+str(i_hPC+1)] = d.obsm['harmpca'][:,i_hPC]
hPC_traits = ['hPC'+str(i_hPC+1) for i_hPC in np.arange(n_hPCs)]
traits = np.concatenate([traits, hPC_traits])

# cell type specific gene expression program
for celltype in celltypes:
    celltype_d = ad.AnnData(d.obsm['harmpca'][d.obs[celltype_col].values==celltype,:20])
    sc.tl.pca(celltype_d)
    for i_celltype_PC in np.arange(n_celltype_PCs):
        d.obs[celltype+"_PC"+str(i_celltype_PC+1)] = np.repeat(np.nan, d.obs.shape[0])
        d.obs.loc[d.obs[celltype_col].values==celltype,celltype+"_PC"+\
                  str(i_celltype_PC+1)]=celltype_d.obsm['X_pca'][:,i_celltype_PC]
celltype_PC_traits = [np.repeat(celltypes,n_celltype_PCs)[i]+"_PC"+\
                      str(np.tile(np.arange(n_celltype_PCs)+1, len(celltypes))[i])\
                      for i in np.arange(n_celltype_PCs*len(celltypes))]
traits = np.concatenate([traits, celltype_PC_traits])

d.obs_to_sample(traits, aggregate = np.mean)

# format for PLINK                                                                                                             
export = d.samplem[traits]
export.insert(0,"#IID", d.samplem.index)
export.to_csv(args.outfile, index=False, sep = "\t", na_rep="NA")
