# Gather expression for selected genes

import numpy as np
import pandas as pd
import cna, os
import scanpy as sc
import argparse
np.random.seed(0)

# Parse Arguments                                                                                                                                          
parser = argparse.ArgumentParser()
parser.add_argument("--lead_snp",type=str)
parser.add_argument("--eQTL_celltype",type=str)
parser.add_argument("--sc_object",type=str)
parser.add_argument("--expr_object",type=str)
parser.add_argument("--csaQTL_celltype",type=str)
parser.add_argument("--out_path",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

celltype=args.eQTL_celltype
csaQTL_celltype=args.csaQTL_celltype

# Import processed data
d = cna.read(args.sc_object)

# Import raw expr
e = sc.read_h5ad(args.expr_object)
e = e[d.obs.index,:] # Subset to QCed cells

# Save raw UMI counts for selected genes 
genes=['TNF', 'IFNG', 'BCL2A1']
genes_ix = [np.where(e.var.index==genes[i])[0][0] for i in np.arange(len(genes))]
genes_exp = pd.DataFrame(e.X[:,genes_ix].todense())
genes_exp.columns=genes
genes_exp.index=e.obs.index
genes_exp.to_csv(args.out_path+csaQTL_celltype+"_"+args.lead_snp+"_csaQTL_test_"+celltype+"_eQTLs_selgene_exp.csv")
