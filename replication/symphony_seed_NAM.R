# Wrapper script to run Symphony and export seed weights per cell for NAM construction

library(argparse)
library(symphony)
library(uwot)
library(harmony)
library(Matrix)
library(RANN)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--celltype",type="character")
parser$add_argument("--cohort",type="character")
parser$add_argument("--ref_folder",type="character")
parser$add_argument("--query_folder",type="character")
parser$add_argument("--out_folder",type="character")
parser$add_argument("--query_batch_vars",type="character") #string list with comma delimiter
args <- parser$parse_args()
cat("\n\n****")
print(args)
cat("****\n\n")

celltype=args$celltype
cohort=args$cohort
ref_folder=args$ref_folder
query_folder=args$query_folder
out_folder=args$out_folder
query_batch_vars = strsplit(args$query_batch_vars,",")[[1]]

# load reference object
reference = readRDS(paste0(ref_folder, celltype, "_symphony_ref.rds"))
ref_gene_ids = read.csv(paste0(ref_folder, celltype,'_vargene_ids.csv')) # Ensemble IDs
reference$vargenes$symbol = ref_gene_ids$gene_id

# load query expression and metadata
query_exp = readMM(paste0(query_folder, cohort, "_", celltype,'_expr.mtx'))
query_genemeta = read.csv(paste0(query_folder,'gene_names_ids.csv'))
query_metadata = read.csv(paste0(query_folder, cohort, "_", celltype,'_cellmeta.csv'))

# reformat query data
query_exp = t(query_exp)
colnames(query_exp) = query_metadata$X # assign cell names (not included in .mtx format)
row.names(query_exp) = query_genemeta$gene_id # label query genes with Ensemble IDs

# run symphony
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_metadata,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 vars = query_batch_vars, # Query batch variables to harmonize over
                 do_normalize = TRUE,  # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP

# export hPC, UMAP embeddings
write.csv(t(query$Z), paste0(out_folder,cohort, "_", celltype,'_hPCs.csv'))
write.csv(query$umap, paste0(out_folder,cohort, "_", celltype,'_umap.csv'))
write.csv(reference$umap, paste0(out_folder,cohort, "_", celltype,'_ref_umap.csv'))

# check whether frequencies of predicted celltypes approximate those in ref dataset
query = knnPredict(query, reference, reference$metadata$celltype)
write.csv(query$meta_data[,c("cell_type_pred_knn","cell_type_pred_knn_prob")], 
          paste0(out_folder,cohort, "_", celltype,'_types.csv'))
write.csv(reference$metadata$celltype, paste0(out_folder,cohort, "_", celltype,'_ref_types.csv'))

# seed query cells into reference NAM using full NN graph
query_nn = nn2(t(reference$Z_corr), t(query$Z), k = 15)
write.csv(query_nn$nn.idx, paste0(out_folder, cohort, "_", celltype,'_in_ref_nngraph_idx.csv'))
write.csv(query_nn$nn.dist, paste0(out_folder,cohort, "_", celltype,'_in_ref_nngraph_dist.csv'))
print("done")



