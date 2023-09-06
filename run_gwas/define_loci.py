# This script defines independent loci and exports the lead SNP for each locus, 
# based on a set of input summary statistics. For all SNPs with associations that pass
# the genome-wide significance threshold, SNPs are sorted in descending order by p-value. 
# The SNP with the lowest p-value is retained as the lead SNP for the first locus. Other 
# SNPs within a 1MB window centered on the lead SNP or with LD>0.8 to the lead SNP 
# (computed using genotyping within the cohort) are removed. Among the remaining SNPs,
# the SNP with the lowest p-value is selected as the lead SNP for the second locus, and so on.
#
# Inputs
#   celltype: title for the GWAS; in our case we distinguish by celltype label
#   outfile: filepath to store GeNA summary statistics file subset to resulting lead SNPs
#   p_thresh: threshold for genome-wide significance
#   gwas_res: filepath to GeNA summary statistics file
#   gtypes: filepath to a paired set of files
#       one with the additional suffix '.vcf.gz' is expected as a gzipped vcf file with genotype information 
#       one with the additional suffix '.DS.vcf.gz' is expected as matching gzipped file with genotype dose only
#   n_gtype_header: number of header lines (to ignore) in the VCF genotype file

import argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments                                                                                                         
parser = argparse.ArgumentParser()
parser.add_argument("--outfile",type=str)
parser.add_argument("--p_thresh",type=float)
parser.add_argument("--gwas_res",type=str)
parser.add_argument("--gtypes",type=str)
parser.add_argument("--celltype",type=str)
parser.add_argument("--n_gtype_header",type=int)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

all_loci = pd.DataFrame({})
window_sz = 1000000 # 1 MB

all_res = pd.read_csv(args.gwas_res, sep = "\t")
all_res['pvalue'] = all_res.P.values
all_res['CHR']=all_res['#CHROM'].values
all_res = all_res.loc[all_res.pvalue<args.p_thresh,:].reset_index(drop=True)

for sel_chr in np.arange(1,23): 
    print(sel_chr)

    # Import GWAS results
    if np.sum(all_res.CHR==sel_chr)==0: 
        print("\tno snps pass genome-wide significance")
        continue
    print("\n")
    res = all_res.loc[all_res.CHR==sel_chr,].reset_index(drop=True)
    res['POS'] = [int(float(res.POS.values[i])) for i in np.arange(res.shape[0])]

    # Import genotype dose
    G = pd.read_csv(args.gtypes+".DS.vcf.gz", sep = "\t", header = None)
    G_meta = pd.read_csv(args.gtypes+".vcf.gz", sep = "\t", header = 1, skiprows = args.n_gtype_header)

    # Pull MAF for reporting
    G_meta['MAF'] = [G_meta.INFO.values[i].split(";")[1].split("=")[1] for i in np.arange(G_meta.shape[0])]
    G_meta.set_index("ID", inplace =True)
    G_chr = [int(G_meta.index[i].split(":")[0]) for i in np.arange(G_meta.shape[0])]
    G = G.loc[G_chr==sel_chr,:]
    G_meta = G_meta.loc[G_chr==sel_chr,:]
    res['MAF'] = G_meta.loc[res.ID.values,'MAF'].values
    res['MAF'] = [float(res.MAF.values[i]) for i in np.arange(res.shape[0])]

    # Sort G and res into decreasing order by p-value
    snp_order = np.argsort(res.pvalue.values)
    res = res.iloc[snp_order,:].reset_index(drop=True)
    G = G.iloc[snp_order,:].reset_index(drop=True)

    # Retain only genome-wide sig SNPs
    keep = np.where(res.pvalue.values<args.p_thresh)[0]
    res = res.iloc[keep,:].reset_index(drop=True)
    G = G.iloc[keep,:].reset_index(drop=True)
    G.set_index([0], inplace = True, drop = True)
    
    while(res.shape[0]>0):
        new = res.iloc[0:1,:]
        new.insert(0, "celltype", args.celltype)
        all_loci = pd.concat([all_loci, new]) # Add lead snp for new locus
        lead_POS = res.POS.values[0]
        res = res.drop(index = res.index[0]) # Remove lead snp
        if res.shape[0]>0: # check if any snps left
            G_corrs = np.abs(G.T.corr()).iloc[1:,:]**2 # Compute R^2 to other remaining snps
            G = G.drop(index = G.index[0]) # Remove lead snp
            i_drop = np.where(G_corrs.iloc[:,0].values>0.8)[0] # Find any other snps in locus
            res = res.drop(index = res.index[i_drop]) # Rm other snps in locus
            G = G.drop(index = G.index[i_drop]) # Rm other snps in locus

            # Also rm anything in 1MB window around lead snp
            in_window = (res.POS > lead_POS-(window_sz/2)) & (res.POS < lead_POS+(window_sz/2))
            i_drop = np.where(in_window)[0]
            res = res.drop(index = res.index[i_drop]) # Rm other snps in locus
            G = G.drop(index = G.index[i_drop]) # Rm other snps in locus
            
all_loci.to_csv(args.outfile, sep = "\t", index = False)
