import pickle, argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--scorename",type=str)
parser.add_argument("--geno_file",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

src_folder= "/data/srlab/lrumker/MCSC_Project/cna-prs/"

# import sumstats
sumstats  = pd.read_csv(src_folder + "sumstats/" + args.scorename+".txt", sep = "\t", header = 14)

# if haplotypes (usually MHC alleles) included in PRS, remove them                                                                                                                                      
if 'is_haplotype' in sumstats.columns:
    print("Removing "+str(np.sum(sumstats.is_haplotype))+" haplotypes; "+str(sumstats.shape[0]-np.sum(sumstats.is_haplotype))+" variants remain")
    sumstats  = sumstats.loc[~sumstats.is_haplotype,:].reset_index(drop = True)

# remove any listed SNPs with effect weight of zero                                                                                                                                                     
sumstats['effect_weight'] = [float(sumstats.effect_weight.values[i]) for i in np.arange(sumstats.shape[0])]
if np.sum(sumstats.effect_weight.values==0)>0:
    print("Removing "+str(np.sum(sumstats.effect_weight.values==0))+" variants with effect weight of 0; "+str(np.sum(sumstats.effect_weight.values!=0))+"variants remain")
    sumstats = sumstats.loc[sumstats.effect_weight!=0,:].reset_index(drop = True)
    
# remove any sex chromosome variants                                                                                                                                                                    
chrs = np.arange(1,23)
chrs = [str(chrs[i]) for i in np.arange(len(chrs))]
autosomal = np.array([str(sumstats.chr_name.values[i]) in chrs for i in np.arange(sumstats.shape[0])])
n_rm = np.sum(~autosomal)
if n_rm > 0:
    sumstats = sumstats.loc[autosomal,:].reset_index(drop = True)
    print("Removing "+str(n_rm)+" non-autosomal variants; "+str(sumstats.shape[0])+" variants remain")

# idenfity variants with available genotypes
sumstats['ID'] = [str(int(sumstats.chr_name[i]))+":"+str(int(sumstats.chr_position[i]))+":"+str(sumstats.other_allele[i])+":"\
                       +str(sumstats.effect_allele[i]) for i in np.arange(sumstats.shape[0])]
sumstats['ID_flip'] = [str(int(sumstats.chr_name[i]))+":"+str(int(sumstats.chr_position[i]))+":"+str(sumstats.effect_allele[i])+":"\
                           +str(sumstats.other_allele[i]) for i in np.arange(sumstats.shape[0])]

# import genotype data
geno = pd.read_csv(args.geno_file, delim_whitespace=True, header = None)
geno['ID'] = geno.iloc[:,0].values

# retain genotyped SNPs
sumstats['found'] = [sumstats.ID.values[i] in geno.ID.values for i in np.arange(sumstats.shape[0])]
sumstats['found_flip'] = [sumstats.ID_flip.values[i] in geno.ID.values for i in np.arange(sumstats.shape[0])]
sumstats['found'] = (sumstats.found.values*1)+(sumstats.found_flip.values*1)
if np.max(sumstats.found.values)>1:
    print("Error: found both a SNP and its flip")

# flip effects as needed
if np.sum(sumstats.found_flip.values)>0:
    print("Switching sign of effect weight (beta) for "+str(np.sum(sumstats.found_flip.values))+" SNPs while flipping their alleles to match available genoypes")
    sumstats.loc[sumstats.found_flip,'effect_weight'] = -1 * sumstats.loc[sumstats.found_flip,'effect_weight'].values
    effect_allele = sumstats.loc[sumstats.found_flip,'effect_allele'].values
    sumstats.loc[sumstats.found_flip,'effect_allele'] = sumstats.loc[sumstats.found_flip,'other_allele'].values
    sumstats.loc[sumstats.found_flip,'other_allele'] = effect_allele

frac_found = np.around(100*(np.sum(sumstats.found)/sumstats.shape[0]),2)
print("Of "+str(sumstats.shape[0])+" SNPs in the published PRS, "+str(np.sum(sumstats.found.values))+" are available ("+str(frac_found)+"%)") 

sumstats['ID'] = [str(int(sumstats.chr_name[i]))+":"+str(int(sumstats.chr_position[i]))+":"+str(sumstats.other_allele[i])+":"\
                       +str(sumstats.effect_allele[i]) for i in np.arange(sumstats.shape[0])]
sumstats.set_index('ID', drop=True, inplace = True)
sumstats = sumstats.loc[geno.ID.values,"effect_weight"]

# save
sumstats.to_csv(src_folder + "sumstats/" + args.scorename+"_refined.txt")

