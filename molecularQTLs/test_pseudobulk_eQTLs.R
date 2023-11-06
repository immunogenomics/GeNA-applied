# Adapted from Kang JB, et al, 2023. Mapping the dynamic genetic regulatory architecture of HLA genes
# at single-cell resolution. medRxiv.
# See: https://github.com/immunogenomics/hla2023/blob/main/scripts/3_pseudobulk/3_run_eQTL_singledataset.R 

library(argparse)
library(data.table)
library(tidyr)
library(dplyr)
library(plyr)
library(Matrix)
library(lme4)
library(Matrix.utils)
library(MASS)
library(glmnet)
require(gdata)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--peer_resid_file",type="character") # path to PEER residuals (.csv)
parser$add_argument("--geno_doses",type="character") # samples X dosages
parser$add_argument("--geno_sample_ids",type="character") # list of sample labels corresponding to geno_doses rows
parser$add_argument("--genes_file",type="character") # path to file with one gene name per row;genes to test
parser$add_argument("--out_path",type="character")
args <- parser$parse_args()
cat("\n\n****")
print(args)
cat("****\n\n")

# load expression data
resid = read.csv(args$peer_resid_file, header = TRUE, sep = ',', row.names = 1) %>% as.matrix() %>% as.data.frame()
test_genes = read.csv(args$genes_file,header = FALSE)[,1]
test_genes = test_genes[test_genes %in% colnames(resid)] # genes that passed QC
resid=resid[,test_genes, drop=FALSE]

# load genotype data
geno <- read.table(args$geno_doses, row.names=1)
geno_ids = read.table(args$geno_sample_ids, header=FALSE)
geno_ids = unname(t(geno_ids)[,1])
colnames(geno) = geno_ids
geno = t(geno[,row.names(resid)])
mod_snp_label<-function(snp_name){
    snp_name=paste0("var",paste(strsplit(snp_name,":")[[1]], collapse="_"), collapse="")
    snp_name=paste(strsplit(snp_name,"<")[[1]], collapse="_")
    return(paste(strsplit(snp_name,">")[[1]], collapse="_"))
}
save_snpnames=colnames(geno)
colnames(geno) = apply(as.matrix(colnames(geno),ncol=1),1,mod_snp_label)

# Add variants to the data frame
data = merge(resid, geno, by = 0)
data$Row.names = NULL # remove extra column added from call to merge

# Run eQTL model
message('Start eQTL modeling')
full_results = NULL
for (gene in test_genes) {
    message('Calculating eQTLs for gene: ', gene)
    results = rbindlist(lapply(colnames(data)[(1+length(test_genes)):ncol(data)], function(variant) { 
        mod = lm(as.formula(paste(gene, ' ~ ', paste(variant, collapse = "+"), sep = "")), data = data)
        
        data.table(variant = variant, 
               gene = gene,
               beta = summary(mod)$coefficients[variant, 'Estimate'], 
               stderr = summary(mod)$coefficients[variant, 'Std. Error'],
               t.val = summary(mod)$coefficients[variant, 't value'],
               p.val = summary(mod)$coefficients[variant, 'Pr(>|t|)'])
    }))
    full_results = rbind(full_results, results)
}

message('Save results')
full_results$variant=rep(save_snpnames, length(test_genes))
write.csv(full_results, paste0(args$out_path, '_pseudobulk_eQTLs.csv'), quote = F)

