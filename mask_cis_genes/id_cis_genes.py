# Identifies the genes within a 2MB window centered on the provided SNP

import argparse
import pandas as pd
import numpy as np
import os
import GTF # from https://gist.github.com/slowkow/8101481
import scanpy as sc
np.random.seed(0)
window_sz = 2000000 # 1 MB upstream, and 1 MB downstream

# Parse Arguments                                                                                                                 
parser = argparse.ArgumentParser()
parser.add_argument("--outfolder",type=str)
parser.add_argument("--celltype",type=str)
parser.add_argument("--gtf_path",type=str)
parser.add_argument("--loci",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

celltype = args.celltype
all_loci = pd.read_table(args.loci)

print("importing GTF")
# gencode.v38lift37.annotation.gtf.gz #gencode.v38.primary_assembly.annotation.gtf
gtf = GTF.dataframe(args.gtf_path)

import pdb;pdb.set_trace()

for i_snp in np.arange(all_loci.shape[0]):

    # chr:pos for SNP of interest
    allele = all_loci.ID.values[i_snp]
    sel_chr = allele.split(":")[0]
    sel_pos = allele.split(":")[1]
    
    gtf_chr = gtf.loc[gtf.seqname.values=="chr"+sel_chr,:].reset_index(drop = True)
    gtf_chr = gtf_chr.loc[gtf_chr.feature=='gene',:].reset_index(drop = True)
    
    # identify genes near snp
    latest_start = int(sel_pos)+int(window_sz/2)
    earliest_end = int(sel_pos)-int(window_sz/2)
    gtf_chr['start'] = [int(gtf_chr.start.values[i]) for i in np.arange(gtf_chr.shape[0])]
    gtf_chr['end'] = [int(gtf_chr.end.values[i]) for i in np.arange(gtf_chr.shape[0])]
    keep = (gtf_chr.end.values>earliest_end) & (gtf_chr.start.values<latest_start)
    gtf_chr = gtf_chr.loc[keep,:].reset_index(drop = True)
    cis_genes = np.unique(gtf_chr.gene_name.values)
    n_cis_genes = len(cis_genes)
    print("Identified "+str(n_cis_genes)+" genes within cis window centered on the variant")
    pd.DataFrame({'cis_genes': cis_genes}).to_csv(args.outfolder+"cis_genes_"+args.celltype+"_"+allele+".tsv", sep = "\t", header = False, index=False)
