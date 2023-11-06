# Adapted from Kang JB, et al, 2023. Mapping the dynamic genetic regulatory architecture of HLA genes
# at single-cell resolution. medRxiv.
# See: https://github.com/immunogenomics/hla2023/blob/main/scripts/3_pseudobulk/2_run_PEER.R

library(argparse)
library(peer)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--expr_file",type="character") # samples X genes, pseudobulked
parser$add_argument("--meta_file",type="character") # samples X covariates
parser$add_argument("--covs_list",type="character") # expected format is cov1,cov2,cov3 all of which must exist as cols in sample_meta
parser$add_argument("--n_peer_factors",type="integer") # number of hidden factors for PEER to use
parser$add_argument("--out_path",type="character")
args <- parser$parse_args()
cat("\n\n****")
print(args)
cat("****\n\n")

# load inputs
expr = read.csv(args$expr_file, header = TRUE, sep = ',', row.names = 1)
meta = read.csv(args$meta_file, header = TRUE, sep = ',', row.names = 1)
covs_list = strsplit(args$covs_list, ",")[[1]]
covs = meta[,covs_list]

message('set up PEER model')
model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model, args$n_peer_factors)
PEER_getNk(model)
PEER_setCovariates(model, as.matrix(covs)) # Add covariates
PEER_setNmax_iterations(model, 10000)

message('Run PEER model')
PEER_update(model)

message('Save results')
residuals = PEER_getResiduals(model)
rownames(residuals) = rownames(expr)
colnames(residuals) = colnames(expr)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)

# Save                        
write.csv(residuals, paste0(args$out_path, '_PEER_residuals.csv'), row.names = TRUE, quote = FALSE)
write.csv(factors, paste0(args$out_path, '_PEER_factors.csv'), row.names = TRUE, quote = FALSE)
write.csv(weights, paste0(args$out_path, '_PEER_weights.csv'), row.names = TRUE, quote = FALSE)
write.csv(precision, paste0(args$out_path, '_PEER_precision.csv'), row.names = TRUE, quote = FALSE)

# Plot convergence
pdf(paste0(args$out_path, "_PEER_plotModel.pdf"), width=8, height=8)
PEER_plotModel(model)
dev.off()

# Plot factor importance
Alpha = PEER_getAlpha(model)
write.csv(Alpha, paste0(args$out_path, '_PEER_alpha.csv'), row.names = TRUE, quote = FALSE)
pdf(paste0(args$out_path, "_PEER_alpha.pdf"), width=8, height=8)
plot(1.0 / Alpha, xlab = "Factors", ylab = "Factor relevance", main="")
dev.off()


