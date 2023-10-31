# Defines cluster-based cell type proportion trait values per individual in the OneK1K dataset  
# For each major cell type (T, B, NK, myeloid, all cells) we compute the fractional abundance
# of each cluster out of 1) all cells, and 2) all cells of the corresponding major cell type.
# We computed pairwise correlations among all these traits and selected one from each pair
# with r-squared > 0.7 for removal in order to maximize the total remaining trait count.
# We also include only traits with nonzero values for at least 200 individuals.
# Values per individual for all retained traits were inverse-normal transformed.

import numpy as np
import pandas as pd
import argparse, cna
import scipy.stats as st

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--outfolder",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

np.random.seed(0)

def inverse_nl_transform(vals):
    N = len(vals)
    obs_vals = pd.DataFrame({"orig_idx": np.arange(N), "val": vals})
    obs_vals = obs_vals.loc[~pd.isna(obs_vals.val),:] # remove nans
    obs_vals = obs_vals.iloc[np.random.choice(np.arange(obs_vals.shape[0]), # Permute order for ties
                                              obs_vals.shape[0], replace = False),:]
    obs_vals = obs_vals.iloc[np.argsort(obs_vals.val),:].reset_index(drop=True) # Sort by values (ties will retain permuted order)
    obs_vals['cumfrac'] = [(obs_vals.index[i]+0.5)/obs_vals.shape[0] for i in np.arange(obs_vals.shape[0])]
    obs_vals['INT_val'] = [st.norm.ppf(obs_vals.cumfrac[i]) for i in np.arange(obs_vals.shape[0])]
    obs_vals.set_index('orig_idx', inplace = True, drop = True)
    INT_vals = pd.DataFrame({"dummy": np.arange(N)}).drop(columns=["dummy"])
    INT_vals = INT_vals.join(obs_vals.loc[:,['INT_val', 'val']])
    return INT_vals

# rm one of each trait pair with R^2>0.7; note that 'NK_frac_allcells' for NK major type already excluded
drop = ['NK_frac_NK', 'CD14Mono_frac_Myeloid', 'Myeloid_frac_allcells', 'CD4 Naive_frac_allcells',
        'CD4 CTL_frac_allcells', 'CD4 TEM_frac_allcells', 'CD8 Naive_frac_allcells',
       'CD8 TCM_frac_allcells', 'CD8 TEM_frac_allcells', 'MAIT_frac_allcells',
       'Treg_frac_allcells', 'gdT_frac_allcells']

# fractional abundance of each minor celltype relative to major celltype
majortypes = ["T", "B", "Myeloid", "NK"]
export = pd.DataFrame()
for denom_celltype in majortypes:
    print(denom_celltype)
    d = cna.read("/data/srlab/lrumker/datasets/onek1k/pheno/"+denom_celltype+".h5ad")
    for celltype in d.obs.celltype.value_counts().index:
        d.obs[celltype+"_frac_"+denom_celltype] = 1*(d.obs.celltype.values==celltype)
        d.obs_to_sample([celltype+"_frac_"+denom_celltype])

    traits = [d.obs.celltype.value_counts().index[i] +"_frac_"+denom_celltype \
              for i in np.arange(d.obs.celltype.value_counts().index.shape[0])]

    # remove traits with values of 0 for 200+ individuals
    traits = np.array(traits)[np.sum(d.samplem[traits]==0,axis=0).values<200] 
    
    # remove one from each trait pair that are strongly correlated
    traits = [traits[i] for i in np.arange(len(traits)) if traits[i] not in drop]

    for sel_trait in traits: # inverse normal transform 
        d.samplem[sel_trait] = inverse_nl_transform(d.samplem[sel_trait].values).INT_val.values

    export = pd.concat([export, d.samplem[traits]], axis=1, join='outer')


denom_celltype = "allcells"
print(denom_celltype)
d = cna.read("/data/srlab/lrumker/datasets/onek1k/pheno/"+denom_celltype+".h5ad")

# fractional abundance of each minor celltype relative to all cells
minortypes = np.unique(d.obs.celltype[d.obs.majortype!='Other']).tolist() # Exclude "other" types
for celltype in minortypes:
    d.obs[celltype+"_frac_"+denom_celltype] = 1*(d.obs.celltype.values==celltype)
    d.obs_to_sample([celltype+"_frac_"+denom_celltype])

# fractional abundance of each major celltype relative to all cells
for celltype in ["T", "B", "Myeloid"]: # NK majortype/all strongly correlates with NK minortype/all, so leave out majortype
    d.obs[celltype+"_frac_"+denom_celltype] = 1*(d.obs.majortype.values==celltype)
    d.obs_to_sample([celltype+"_frac_"+denom_celltype])
    
traits_minor = [minortypes[i] +"_frac_"+denom_celltype for i in np.arange(len(minortypes))]
traits_major = ['T_frac_allcells', 'B_frac_allcells', 'Myeloid_frac_allcells']
traits = np.concatenate([traits_minor, traits_major])

# remove traits with values of 0 for 200+ individuals                                                                                
traits = np.array(traits)[np.sum(d.samplem[traits]==0,axis=0).values<200]

# remove one from each trait pair that are strongly correlated
traits = [traits[i] for i in np.arange(len(traits)) if traits[i] not in drop]

for sel_trait in traits: # inverse normal transform                                                                                    
    d.samplem[sel_trait] = inverse_nl_transform(d.samplem[sel_trait].values).INT_val.values
export = pd.concat([export, d.samplem[traits]], axis=1, join='outer')

# format for PLINK
fid = [export.index[i].split("_")[0] for i in np.arange(export.shape[0])]
iid = [export.index[i].split("_")[1] for i in np.arange(export.shape[0])]
export.insert(0,"IID", iid)
export.insert(0,"FID", fid)
export.to_csv(args.outfolder+"cluster_traits.tsv", index=False, sep = "\t", na_rep="NA")

# export covs for PLINK                                                                                                      
export_covs = d.samplem.loc[:,["sex_M", "age", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5", "gPC6"]]
fid = [d.samplem.index[i].split("_")[0] for i in np.arange(d.samplem.shape[0])]
iid = [d.samplem.index[i].split("_")[1] for i in np.arange(d.samplem.shape[0])]
export_covs.insert(0,"IID", iid)
export_covs.insert(0,"FID", fid)
export_covs.to_csv(args.outfolder+'covs.tsv', index=False, sep = "\t", na_rep = "NA")
