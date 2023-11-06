import pickle, argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--scorename",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

src_folder= "/data/srlab/lrumker/MCSC_Project/cna-prs/sumstats/"
try_flips = True

sumstats  = pd.read_csv(src_folder +args.scorename+".txt", sep = "\t", header = 14)

print(str(sumstats.shape[0])+ " variants in PRS")

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

# in geno data, snps labeled as chr:pos:ref:alt and record dose of alt allele
# so we take "effect allele" in sumstats to be alt allele
sumstats['SNP_ID'] = [str(int(sumstats.chr_name[i]))+":"+str(int(sumstats.chr_position[i]))+":"+str(sumstats.other_allele[i])+":"+str(sumstats.effect_allele[i]) for i in np.arange(sumstats.shape[0])]
if try_flips: # consider flipped snps
    sumstats['SNP_FLIP'] = [str(int(sumstats.chr_name[i]))+":"+str(int(sumstats.chr_position[i]))+":"+str(sumstats.effect_allele[i])+":"+str(sumstats.other_allele[i]) for i in np.arange(sumstats.shape[0])]
    pd.DataFrame({"SNP":np.concatenate((sumstats.SNP_ID.values, sumstats.SNP_FLIP.values))}).to_csv(src_folder+"/snplists/"+args.scorename+"_snplist.txt", header = False, index = False)
else:
    pd.DataFrame({"SNP":sumstats.SNP_ID.values}).to_csv(src_folder+"/snplists/"+args.scorename+"_snplist.txt", header = False, index = False)
    
