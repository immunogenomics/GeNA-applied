# Adapted from Kang JB, et al, 2023. Mapping the dynamic genetic regulatory architecture of HLA genes
# at single-cell resolution. medRxiv.
# See https://github.com/immunogenomics/hla2023/blob/main/scripts/3_pseudobulk/1_pseudobulk_normalize.R

library(argparse)
library(tidyr)
library(dplyr)
library(plyr)
library(Matrix)
library(Matrix.utils)
library(singlecellmethods)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--celltype_expr",type="character")
parser$add_argument("--celltype_meta",type="character")
parser$add_argument("--genenames",type="character")
parser$add_argument("--out_path",type="character")
parser$add_argument("--sample_meta",type="character")
args <- parser$parse_args()
cat("\n\n****")
print(args)
cat("****\n\n")

# load expression and metadata
exprs_raw = readMM(args$celltype_expr) #raw UMI counts, mm file format
exprs_genes = read.csv(args$genenames, header = FALSE)$V1
cell_meta = read.csv(args$celltype_meta, row.names=1)
row.names(exprs_raw) = row.names(cell_meta)
colnames(exprs_raw) = exprs_genes
exprs_raw = t(exprs_raw) # result is genes X cells dimensions

# consider only samples with at least 25 cells
cc_df = data.frame(table(cell_meta$individual))
colnames(cc_df) = c("ID", "cell_count")
cc_df = cc_df[cc_df$cell_count>24,]
keep_sample <-function(sel_ID){
   return(sel_ID %in% cc_df$ID)
}
cell_meta$keep = apply(as.matrix(cell_meta$individual, ncol=1), 1, keep_sample)
exprs_raw = exprs_raw[,cell_meta$keep]
cell_meta = cell_meta[cell_meta$keep,]

# Make pseudobulk profiles by single-cell log-normalization then mean aggregation across cells
exp_norm = normalizeData(exprs_raw, method = "log")
exp_norm = exp_norm[which(rowSums(exp_norm) != 0), ] # remove genes that are 0 across all cells
pseudobulk_scnorm_exp_sum = aggregate.Matrix(t(exp_norm), as.factor(cell_meta$individual), fun = 'sum') # samples x genes
pseudobulk_scnorm_exp_mean = pseudobulk_scnorm_exp_sum / count(cell_meta$individual)$freq

# Include only genes with non-zero expression in > half of the samples
pseudobulk_scnorm_exp_mean = pseudobulk_scnorm_exp_mean[, colSums(pseudobulk_scnorm_exp_mean > 0) > .5 * nrow(pseudobulk_scnorm_exp_mean)]

# Subset sample metadata to only the included samples
sample_meta=read.csv(args$sample_meta, row.names=1)
sample_meta=sample_meta[,c("gPC1","gPC2","gPC3","gPC4","gPC5","gPC6")]
cell_meta_bySample = cell_meta[!duplicated(cell_meta$individual),]
row.names(cell_meta_bySample) = cell_meta_bySample$individual
cell_meta_bySample = cell_meta_bySample[row.names(sample_meta),c("age", "sex")]
sample_meta['age'] = cell_meta_bySample$age
sample_meta['sex'] = cell_meta_bySample$sex-1
sample_meta_subset = sample_meta[rownames(pseudobulk_scnorm_exp_mean), ]
write.csv(sample_meta_subset, paste0(args$out_path, '_samples_meta.csv'), row.names = TRUE, quote = FALSE) # Write sample meta
if(all(rownames(pseudobulk_scnorm_exp_mean) != rownames(sample_meta_subset))) {
    message('Error: ordering of samples inconsistent.')}

# Write normalized result
saveRDS(pseudobulk_scnorm_exp_mean, paste0(args$out_path, '_samplesXgenes_norm.rds'))

# Rank-based inverse normal transformation
exp_int = apply(pseudobulk_scnorm_exp_mean, 2, function(x) {qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )})

# Write inverse normal transformed result
write.csv(exp_int, paste0(args$out_path, '_samplesXgenes_norm_invnt.csv'), row.names = TRUE, quote = FALSE)