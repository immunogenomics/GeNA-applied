import numpy as np
import pandas as pd
import cna, os
import scanpy as sc
import argparse
np.random.seed(0)

# Parse Arguments                                                                                                                                          
parser = argparse.ArgumentParser()
parser.add_argument("--eQTL_celltype",type=str)
parser.add_argument("--sc_object",type=str)
parser.add_argument("--expr_object",type=str)
parser.add_argument("--pbulk_res_folder",type=str)
parser.add_argument("--p_thresh",type=float)
parser.add_argument("--out_path",type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

celltype=args.eQTL_celltype

# Import processed data
d = cna.read(args.sc_object)

# Import raw expr
e = sc.read_h5ad(args.expr_object)
e = e[d.obs.index,:] # Subset to QCed cells

# Save gene expression PCs
ePCs = pd.DataFrame(d.obsm['X_pca'][:,:5])
ePCs.index=d.obs.index
ePCs.columns = ["ePC"+str(i+1) for i in np.arange(5)]
ePCs.to_csv(args.out_path+celltype+"_ePCs.csv")

# Save cell metadata
cell_meta = d.obs.loc[:,['id', 'batch', 'sex', 'age', 'nCount_RNA', 'percent.mt', 'celltype']]
for gPC in ['gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5', 'gPC6']:
   cell_meta[gPC] = d.samplem.loc[d.obs.id.values, gPC].values
cell_meta.to_csv(args.out_path+celltype+"_cellmeta.csv")

# Identify candidate eGenes based on pseudobulk results
files=os.listdir(args.pbulk_res_folder)
files = [files[i] for i in np.arange(len(files)) if "csaQTL_test_"+celltype+"_eQTLs_pseudobulk_eQTLs" in files[i]]

for res_file in files:
   res = pd.read_csv(args.pbulk_res_folder+res_file, index_col=0)
   csaQTL_celltype=res_file.split("_")[0]
   lead_snp=res_file.split("_")[1]
   genes = res.loc[np.logical_and(res.variant==lead_snp, res['p.val']<args.p_thresh),'gene'].values

   if len(genes)>0:
      # Save raw UMI counts for selected genes 
      genes_ix = [np.where(e.var.index==genes[i])[0][0] for i in np.arange(len(genes))]
      genes_exp = pd.DataFrame(e.X[:,genes_ix].todense())
      genes_exp.columns=genes
      genes_exp.index=e.obs.index
      genes_exp.to_csv(args.out_path+csaQTL_celltype+"_"+lead_snp+"_csaQTL_test_"+celltype+"_eQTLs_selgene_exp.csv")
