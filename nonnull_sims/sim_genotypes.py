
# This script generates simulated genotypes with the specified MAF and                                                                  
# true associations to the input cell state abundance traits.                                                                           
# To introduce noise and vary the amount of variance explained in the trait                                                             
# by the simulated genotype, for each phenotype we permuted genotype values for                                                         
# 1%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80% and 100% of samples. For each                                                              
# of these ten tiers of sample count to permute, the script will generate 'nsim'                                                        
# genotype permutations (simulates). For each simulate, we selected the samples of 
# desired count to permute at random with equal probability among all samples. We
# permuted genotype values among the selected samples at random.
# As the count of samples included in the permutation increases, 
# the phenotypic variance explained by the resulting simulated SNP decays. 

import numpy as np
import pandas as pd
import argparse
np.random.seed(0)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--celltype",type=str)
parser.add_argument("--outfile_genos",type=str)
parser.add_argument("--outfile_meta",type=str)
parser.add_argument("--trait_file",type=str)
parser.add_argument("--nsim",type=int)
parser.add_argument("--maf",type=float)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

maf = args.maf
celltype = args.celltype
n_sim = args.nsim

# load traits
traits = pd.read_table(args.trait_file)
print("Simulating genotypes corresponding to the following traits:")
print(traits.columns[1:])
ids = traits.loc[:,['#IID']]
traits.drop(columns=['#IID'], inplace = True)
traits.columns = ["T"+str(i+1) for i in np.arange(traits.shape[1])]

# objects to store simulated genotypes
gtype_meta_all=pd.DataFrame({})
sim_gtypes_all=pd.DataFrame({})

# genotype frequencies for this maf at HW equilibrium
q = maf
p = 1 - q
freq_AA = q**2
freq_RR = p**2

# simulations per trait
for trait in traits.columns:
 
    sim_gtypes=np.array([[]])

    # trait values in rank order (save original order)
    trait_vals = np.array(traits[trait].values)
    trait_order = np.arange(len(trait_vals))
    rank_order = np.argsort(-trait_vals)
    trait_vals = trait_vals[rank_order]
    trait_order = trait_order[rank_order]
    trait_isna = pd.isna(trait_vals)
    trait_vals = trait_vals[~trait_isna]
    trait_order = trait_order[~trait_isna]

    # genotype values at this MAF for maximal variance explained
    N=len(trait_vals)
    N_RR = int(np.around(freq_RR*N,0))
    N_AA = int(np.around(freq_AA*N,0))
    N_RA = int((N-N_AA)-N_RR)
    gtypes = np.concatenate([np.concatenate([np.repeat(0, N_RR), np.repeat(1, N_RA)]), np.repeat(2, N_AA)])

    noise_levels = [int(len(trait_vals)*i) for i in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]]
    
    # iterate through noise levels
    for n_shuffle in noise_levels:    

        # make simulated genotypes at the given noise level
        for n_iter in np.arange(n_sim):
            i_to_shuffle = np.random.choice(np.arange(len(trait_vals)), size = n_shuffle, replace = False)
            i_post_shuffle = np.random.choice(i_to_shuffle, size = n_shuffle, replace = False)
            i_gtypes = np.arange(len(gtypes))
            i_gtypes[i_post_shuffle] = i_to_shuffle
            gtypes_shuffled = gtypes[i_gtypes]

            new_gtypes = np.repeat(np.nan, traits.shape[0])
            new_gtypes[trait_order] = gtypes_shuffled

            if sim_gtypes.shape[1]==0:
                sim_gtypes=new_gtypes.reshape(-1,1)
            else:
                sim_gtypes = np.hstack([sim_gtypes, new_gtypes.reshape(-1,1)])

    # compute variance explained for each simulated genotype
    sim_gtypes_mod = pd.DataFrame(sim_gtypes)
    sim_gtypes_mod.columns=['G'+str(i+1) for i in np.arange(sim_gtypes.shape[1])]
    sim_gtypes_mod[trait] = traits[trait].values

    # store simulate meta information
    gtype_meta = pd.DataFrame({"varexpl": sim_gtypes_mod.corr().loc[trait][:-1]**2,
                              'trait': np.repeat(trait, sim_gtypes.shape[1]),
                              'n_shuffle': np.repeat(noise_levels, n_sim)})

    gtype_meta_all = pd.concat([gtype_meta_all, gtype_meta])
    sim_gtypes_all = pd.concat([sim_gtypes_all, pd.DataFrame(sim_gtypes)], axis=1)

# Save simulated genotypes as vcf
tot_simulates = sim_gtypes_all.T.shape[0]
gtype_meta_all['n_simulate'] = np.tile(np.arange(n_sim)+1,len(noise_levels)*traits.shape[1])
gtype_meta_all['sim_id'] = [gtype_meta_all.trait.values[i]+"_shuffle"+str(gtype_meta_all.n_shuffle.values[i])+\
                            "_sim"+str(gtype_meta_all.n_simulate.values[i]) for i in np.arange(gtype_meta_all.shape[0])]
gtype_meta_all.to_csv(args.outfile_meta, header = True, index = False, sep = "\t")

dummy_metadata = pd.DataFrame({'#CHROM': np.repeat(1, tot_simulates),
                          'POS': np.arange(tot_simulates)+100,                               
                          'ID': gtype_meta_all['sim_id'].values,
                          'REF': np.repeat('A', tot_simulates),
                          'ALT': np.repeat('T', tot_simulates),
                          'QUAL': np.repeat('.', tot_simulates),
                          'FILTER': np.repeat('PASS', tot_simulates),
                          'INFO': np.repeat('AF=0.35081;MAF=0.35081;R2=0.82863', tot_simulates),
                          'FORMAT': np.repeat('GT:DS', tot_simulates)})

# Update genotypes to GT:DS format
sim_gtypes_all = sim_gtypes_all.replace(0, '0|0:0')
sim_gtypes_all = sim_gtypes_all.replace(1, '1|0:1')
sim_gtypes_all = sim_gtypes_all.replace(2, '1|1:2')
sim_gtypes_all = sim_gtypes_all.T

# add real sample labels to simulated genotypes
sim_gtypes_all.columns = ids['#IID'].values

# make VCF body
sim_gtypes_all.reset_index(inplace = True, drop = True)
dummy_metadata.reset_index(inplace = True, drop = True)
export = pd.concat((dummy_metadata,sim_gtypes_all),axis=1)
export.to_csv(args.outfile_genos, header = True, index = False, sep = "\t", na_rep = ".")
