# Generates a single-cell object that drops information for any variable genes from the discovery data
# object that lie within a 2MB cis window around the lead SNP 
# Defines the phenotype most correlated with the lead SNP in the masked object

import argparse
import pandas as pd
import numpy as np
import cna, os, pickle
import scanpy as sc
import anndata as ad
import multianndata as mad
import harmonypy as hm
covs = ['age','sex_M', 'gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5', 'gPC6']
src_folder = "/data/srlab/lrumker/datasets/onek1k/pheno/"
np.random.seed(0)
window_sz = 2000000 # 1 MB upstream, and 1 MB downstream

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--cis_genes",type=str)
parser.add_argument("--celltype",type=str)
parser.add_argument("--lead_snp",type=str)
parser.add_argument("--gtypes",type=str)
parser.add_argument("--gtype_samples",type=str)
parser.add_argument("--cna_res",type=str)
parser.add_argument("--out_sc_object",type=str)
parser.add_argument("--sel_k",type=int)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

celltype = args.celltype
allele = args.lead_snp
sel_chr = int(allele.split(":")[0])
cis_genes = pd.read_csv(args.cis_genes, header = None).iloc[:,0].values

print("Reading in data for celltype: "+celltype)
d = cna.read('/data/srlab/lrumker/datasets/onek1k/pheno/'+celltype+'.h5ad')
d.scale_variance = False
d.use_R2 = False
d.count_factor = 0

print("Reading in data for chromosome: "+str(sel_chr))
G = pd.read_csv(args.gtypes, sep = "\t", index_col = 0, header = None)
G_samples = pd.read_table(args.gtype_samples, header = None, sep = "\t")
G.columns = G_samples.iloc[:,0].values
G.index = [G.index[i].strip() for i in np.arange(G.shape[0])]
d.samplem = d.samplem.join(G.T)

found = [cis_genes[i] in d.var.index for i in np.arange(len(cis_genes))]
print(str(len(cis_genes))+" genes found in a 2MB cis window centered on the SNP.")
print(str(np.sum(found))+" of these genes were included in original variable genes set, now removed:")
print(cis_genes[found])
keep_vargenes = [i for i in np.arange(d.var.shape[0]) if d.var.index[i] not in cis_genes]
genenames = d.var.index[keep_vargenes]
d_mc = ad.AnnData(X = d.X[:,keep_vargenes], obs = d.obs) # expr already logarithmized, normalized, scaled
d_mc.var.index = genenames

sc.tl.pca(d_mc, svd_solver='arpack') # PCA
if celltype == "allcells":
    ho = hm.run_harmony(d_mc.obsm['X_pca'][:,:20], d_mc.obs, ['pool'], max_iter_harmony = 50, theta = 2)
else:
    ho = hm.run_harmony(d_mc.obsm['X_pca'][:,:20], d_mc.obs, ['pool'], nclust = 50, sigma = 0.2, max_iter_harmony = 50, theta = 2)
d_mc.obsm['harmpca'] = ho.Z_corr.T

print("nn graph")
sc.pp.neighbors(d_mc, use_rep = 'harmpca') # nn graph                                                                            

d_mc = mad.MultiAnnData(d_mc, sampleid='id')
d_mc.use_R2 = False
d_mc.scale_variance = False
d_mc.count_factor = 0

d_mc.obs_to_sample(['sex', 'age', 'batch'])
d_mc.samplem['sex_M'] = (d_mc.samplem.sex==1)*1 # From 1 vs 2 to boolean
d_mc.samplem = d_mc.samplem.drop(columns = ['sex'])
meta = pd.read_csv(src_folder+"sample_meta.csv", index_col = 0)
d_mc.samplem = d_mc.samplem.join(meta.loc[:,['gPC1', 'gPC2', 'gPC3', 'gPC4', 'gPC5', 'gPC6']])
d_mc.samplem = d_mc.samplem.join(G.T)

cna.tl.nam(d_mc, batches=d_mc.samplem.batch, covs=d_mc.samplem[covs])
d_mc.write(args.out_sc_object)

res = cna.tl.association(d_mc, d_mc.samplem[allele], covs = d_mc.samplem[covs], batches = d_mc.samplem.batch,
                         ks = [args.sel_k])
print("vargene correlations")
vargene_cors = []
for i_gene in np.arange(d_mc.var.shape[0]):
    vargene_cors.append(np.corrcoef(d_mc.X[res.kept, i_gene], res.ncorrs)[0,1])
res.vargene_cors=pd.DataFrame({'gene':d_mc.var.index, 'cor': vargene_cors})

pickle.dump(res, open(args.cna_res, 'wb'))

