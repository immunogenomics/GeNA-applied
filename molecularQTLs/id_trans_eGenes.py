# Identifies the trans-genes (outside 2MB window centered on the provided SNP) with even weakly
# suggestive (p<1e-4) eQTL associations to the provided csaQTL lead SNP in a pseudobulk model

import argparse
import pandas as pd
import numpy as np
np.random.seed(0)

# Parse Arguments                                                                                                                 
parser = argparse.ArgumentParser()
parser.add_argument("--outfile",type=str)
parser.add_argument("--celltype",type=str)
parser.add_argument("--lead_snp",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

# Note that not all vargenes passed PEER QC for pseudobulk eQTL testing
res = pd.read_csv("/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/trans_eQTL_candidates/pseudobulk/"+
              args.celltype+"_"+args.lead_snp+"_csaQTL_test_eQTLs_pseudobulk_eQTLs.csv", index_col = 0)
res = res.loc[res.variant.values==args.lead_snp,:].reset_index(drop=True)
res.set_index('gene', inplace = True)

# Suggestive eQTL p-value threshold
res = res.loc[res['p.val']<1e-4,:]

# If any candidates, exclude cis genes
if res.shape[0]==0:
    pd.DataFrame({"eGene": []}).to_csv(args.outfile, header=False, index = False)
else:
    cis_genes = pd.read_csv("/data/srlab/lrumker/MCSC_Project/cna-qtl/mask_cis/results/cis_genes_"+\
                    args.celltype+"_"+args.lead_snp+".tsv", header = None).iloc[:,0].values
    is_cis = [res.index[i] in cis_genes for i in np.arange(res.shape[0])]
    res = res.loc[~np.array(is_cis),:]
    pd.DataFrame({"eGene": res.index}).to_csv(args.outfile, header=False, index = False)
