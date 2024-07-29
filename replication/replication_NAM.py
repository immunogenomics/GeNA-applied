import argparse
import pandas as pd
import numpy as np
import cna
import scanpy as sc
import pp, os
import statsmodels.api as sm
from projection import *
np.random.seed(0)

# Parse arguments                                                                                                                       
parser = argparse.ArgumentParser()
parser.add_argument("--celltype",type=str)
parser.add_argument("--cohort",type=str)
parser.add_argument("--save_NAM",type=bool,default=False)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

ref_folder = "/data/srlab/lrumker/datasets/onek1k/pheno/"
src_folder = "/data/srlab/lrumker/datasets/perez_sle/pheno/"
minor_types = {"NK": ["NK"], "T": ['T4',  'T8'], "Myeloid": ["cM", "ncM", 'cDC', 'pDC'], "B": ['B',  'PB']}

celltype=args.celltype
cohort=args.cohort
covs_perez = ['age', 'sle']

# Load discovery dataset
r = cna.read(ref_folder+celltype+".h5ad")
covs = ['age']
attr = 'sex_M'
cna.tl.nam(r, batches=r.samplem['batch'], covs = r.samplem[covs], ks=[r.samplem.shape[0]]) #to store all nampcs
res = cna.tl.association(r, r.samplem[attr], batches=r.samplem['batch'], covs = r.samplem[covs])

d = cna.read(src_folder+"symphony_output/projection_examples/"+cohort+"_"+celltype+".h5ad")

# Import NAM seed information
pred_idx = (pd.read_csv(src_folder+'symphony_output/projection_examples/'+cohort+"_"+celltype+"_in_ref_nngraph_idx.csv", 
                        index_col = 0).reset_index(drop = True)-1) # R to python indexing
pred_dist = pd.read_csv(src_folder+'symphony_output/projection_examples/'+cohort+"_"+celltype+"_in_ref_nngraph_dist.csv", 
                        index_col = 0).reset_index(drop = True)
pred_sim = 1/pred_dist # similarity

# Construct the NAM
print("seeding NAM")
NAM = seed_nam(r, d, pred_idx, pred_sim, sampleid="id")
print("diffusing NAM")
NAM = nam_diffuse_from_seed(r, NAM.T)
if args.save_NAM:
    pd.DataFrame(NAM).to_csv("/data/srlab/lrumker/MCSC_Project/cna-qtl/gwas_replication/replication_sex_"+celltype+"_"+cohort+"_NAM.csv")
        
# Project phenotype to replication dataset
k = res.k
D = np.identity(len(r.uns['NAM_svs'][:k]))*r.uns['NAM_svs'][:k]
D_I = np.linalg.inv(D)
mask = np.abs(res.ncorrs)>res.fdr_10p_t
sampXpc = NAM[:,r.uns['keptcells']][:,mask].dot(r.uns['NAM_nbhdXpc'].iloc[:,:k].loc[mask,:]).dot(D_I)
est_pheno = np.dot(sampXpc[:,:k], res.beta).reshape(-1,)
features = sm.add_constant(d.samplem[covs_perez+[attr]])
features.loc[:,attr] = features.loc[:,attr].values.astype(float)
linmod = sm.OLS(est_pheno.astype(float), features).fit() 

# P-value from a one-tailed test for concordant-direction association
pval = linmod.pvalues.loc[attr]
beta = linmod.params[attr]
if beta > 0:
    pval = pval/2 
else: 
    pval = 1-(pval/2)

new = pd.DataFrame({"cohort":[cohort], "celltype":[celltype], "attr":[attr], 'P':[pval],
                    'beta':[linmod.params[attr]], 'stderr':[linmod.bse[attr]]})

new.to_csv("/data/srlab/lrumker/MCSC_Project/cna-qtl/gwas_replication/replication_sex_"+celltype+"_"+cohort+".csv") # save 
