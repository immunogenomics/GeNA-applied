import pickle, argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# MHC coordinates for hg19
HLA_START=28477797
HLA_STOP=33448354

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--scorename",type=str)
parser.add_argument("--remove_MHC",type=bool, default=False)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

src_folder= "/data/srlab/lrumker/MCSC_Project/cna-prs/"

# import data objects
sumstats  = pd.read_csv(src_folder +"sumstats/" + args.scorename+"_refined.txt")
geno = pd.read_csv(src_folder + "results/geno/" + args.scorename+".all_chr.DS.vcf.gz", sep = "\t", index_col = 0, header = None)
geno.index = [geno.index.values[i].strip() for i in np.arange(geno.shape[0])]
samples = pd.read_csv(src_folder + "results/geno/sample_list.txt", header = None)
geno.columns = samples.iloc[:,0]

# Ignore MHC
if args.remove_MHC:
    sumstats['POS'] = [int(sumstats.ID.values[i].split(":")[1]) for i in np.arange(sumstats.shape[0])]
    sumstats['CHR'] = [int(sumstats.ID.values[i].split(":")[0]) for i in np.arange(sumstats.shape[0])]
    keepsnp = (sumstats.POS < HLA_START) | (sumstats.POS > HLA_STOP) | (sumstats.CHR !=6 )
    print("Removing "+str(np.sum(~keepsnp))+" MHC SNPs; "+str(np.sum(keepsnp))+" variants remain")
    sumstats = sumstats.loc[keepsnp,:].reset_index(drop=True)
    geno = geno.loc[sumstats.ID.values,:]

PRS = np.dot(geno.T, sumstats.effect_weight.values.reshape(-1,1))

# Save 
filename = src_folder+"results/PRS/"+args.scorename+"_PRS.csv"
if args.remove_MHC:
    filename = src_folder+"results/PRS/"+args.scorename+"_PRS_noMHC.csv"
pd.DataFrame({"PRS": PRS.reshape(-1,) , "IDs": geno.columns}).to_csv(filename)

