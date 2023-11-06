import numpy as np 
import pandas as pd
import scanpy as sc
import cna, os
import anndata as ad
import multianndata as mad
import scipy.stats as st
import gc

# Seeds density of cells from dataset d in neighborhoods
# for dataset r using embedding of cells from d in nn graph
# for r, captured by idx (indices of direct neighbors) and sim
# (similarities of the cell to each of the direct neighbors)
def seed_nam(r, d, idx, sim, sampleid="id", verbose = True):
    NAM = np.zeros((0,r.shape[0]))
    for j, ind_cov in enumerate(d.samplem.index):
        if verbose:
            print("Projecting sample "+str(ind_cov)+" ("+str(j+1)+" of "+str(d.samplem.shape[0])+")")
        A = np.zeros((np.sum(d.obs[sampleid].values==ind_cov), r.shape[0]))
        sample_cells = idx.loc[d.obs[sampleid].values==ind_cov,:]
        sample_sim = sim.loc[d.obs[sampleid].values==ind_cov,:]
        sample_cells.reset_index(inplace = True, drop = True)
        sample_sim.reset_index(inplace = True, drop = True)
        for i, row in enumerate(sample_cells.iterrows()):
            A[row] = sample_sim.iloc[i,:].values
        A = A / A.sum(axis=1)[:,None] # each cell sums to 1
        NAM = np.vstack([NAM, A.sum(axis=0)])
    NAM = NAM/NAM.sum(axis=1)[:,None] # each sample sums to 1
    return NAM

# creates a neighborhood abundance matrix
# uses neighborhood connectivities given in ref_data
# takes an already-seeded NAM ("S")
def nam_diffuse_from_seed(ref_data, S, nsteps=None, maxnsteps=15):
    def R(A, B):
        return ((A - A.mean(axis=0))*(B - B.mean(axis=0))).mean(axis=0) \
            / A.std(axis=0) / B.std(axis=0)
    
    C = S.sum(axis=0)
    
    prevmedkurt = np.inf
    old_s = np.zeros(S.shape)
    for i, s in enumerate(cna.tl.diffuse_stepwise(ref_data, S, maxnsteps=maxnsteps)):
        medkurt = np.median(st.kurtosis(s/C, axis=1))
        print('\tmedian kurtosis:', medkurt+3)
        if np.sum(s.sum(axis=1)==0)==0 and i>0:
            R2 = R(s, old_s)**2
            print('\t20th percentile R2(t,t-1):', np.percentile(R2, 20))
        old_s = s
        if nsteps is None:
            if prevmedkurt - medkurt < 3 and i+1 >= 3:
                print('stopping after', i+1, 'steps')
                break
            prevmedkurt = medkurt
        elif i+1 == nsteps:
            break
        gc.collect()
    snorm = (s / C).T
    return snorm
