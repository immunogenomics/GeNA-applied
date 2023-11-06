# Adapted from single-cell Poisson mixed effects eQTL model without cell state interaction from Nathan A, et al, 2023. Single-cell eQTL models reveal dynamic T cell state dependence of disease loci. Nature.
# Specifically, we apply the single-cell linear mixed effects eQTL model without an eQTL <-> cell state interaction
# See: https://github.com/immunogenomics/sceQTL/blob/main/scripts/singlecell/linear_nostate.R

library(argparse)
library(lme4)
library(Matrix)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--lead_snp",type="character")
parser$add_argument("--src_filepath",type="character")
parser$add_argument("--res_filepath",type="character")
parser$add_argument("--celltype_expr",type="character") #candidate eQTL celltype
parser$add_argument("--celltype_geno",type="character") #csaQTL celltype
parser$add_argument("--geno",type="character")
parser$add_argument("--geno_ids",type="character")
parser$add_argument("--p_thresh",type="double")
args <- parser$parse_args()
cat("\n\n****")
print(args)
cat("****\n\n")

# check if any candidate eGenes
expr_file=paste0(args$src_filepath,args$celltype_geno,'_',args$lead_snp,'_csaQTL_test_',args$celltype_expr,'_eQTLs_selgene_exp.csv')
if(!file.exists(expr_file)){
    print(paste0("No candidate eGenes for ",expr_file))
    quit()
}

# load phenotype and covariate data
exprs_raw = read.csv(expr_file, row.names = 1) #raw UMI counts
pca_res = read.csv(paste0(args$src_filepath, args$celltype_expr, "_ePCs.csv"), row.names = 1) # gene expression PCs
cell_meta = read.csv(paste0(args$src_filepath,args$celltype_expr, "_cellmeta.csv")) # cell and donor covariates

# test per eGene candidate
for(gene in colnames(exprs_raw)){
    print(gene)
    data = cbind(cbind(exprs_raw, pca_res), cell_meta)
    data[,gene] = as.numeric(data[,gene]) 
    data$id = factor(data$id)
    data$age = scale(data$age)
    data$nCount_RNA = scale(log(data$nCount_RNA)) # nUMI

    # load genotype data
    geno <- read.table(args$geno, row.names=1)
    geno_ids = read.table(args$geno_ids)
    geno_ids = as.character(geno_ids[1,])
    colnames(geno) = geno_ids
    geno = geno[,colnames(geno) %in% unique(data$id)] # only donors that passed QC

    data['E'] = data[,gene]

    # Just test lead snp, stop if assoc below threshold
    G_snp = data.frame("G" = as.numeric(as.character(geno[rownames(geno)==args$lead_snp, match(data$id, colnames(geno))])))
    mod_data = cbind(data, G_snp)
    tryCatch({
        full_model <- lme4::glmer(formula = E~G+(1|id)+(1|batch)+age+sex+nCount_RNA+percent.mt+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+ePC1+ePC2+ePC3+ePC4+ePC5,
                              family = "poisson", nAGQ=0, data= mod_data, control = glmerControl(optimizer = "nloptwrap"))
        null_model <- lme4::glmer(formula = E~(1|id)+(1|batch)+age+sex+nCount_RNA+percent.mt+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+ePC1+ePC2+ePC3+ePC4+ePC5,
                                  family = "poisson", nAGQ=0, data= mod_data, control = glmerControl(optimizer = "nloptwrap"))
        model_lrt <- anova(null_model, full_model)
        res = data.frame("SNP" = args$lead_snp, "GENE" = gene,
                 "BETA" = summary(full_model)$coefficients[2,][1], #G beta
                "SE" = summary(full_model)$coefficients[2,][2], #G se
                "P" = model_lrt$`Pr(>Chisq)`[2])
        }, error=function(cond){return(NA)})

    if(res$P>args$p_thresh){
	print(paste0("Association for lead SNP has p=",res$P," so the remaining loci will not be tested."))
	next
    }

    all_res = data.frame()
    for( snp in rownames(geno)){
    	print(which(rownames(geno)==snp)/nrow(geno))
    	G_snp = data.frame("G" = as.numeric(as.character(geno[rownames(geno)==snp, match(data$id, colnames(geno))])))
    	mod_data = cbind(data, G_snp)
    	tryCatch({
	    full_model <- lme4::glmer(formula = E~G+(1|id)+(1|batch)+age+sex+nCount_RNA+percent.mt+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+ePC1+ePC2+ePC3+ePC4+ePC5,
                              family = "poisson", nAGQ=0, data= mod_data, control = glmerControl(optimizer = "nloptwrap"))
            null_model <- lme4::glmer(formula = E~(1|id)+(1|batch)+age+sex+nCount_RNA+percent.mt+gPC1+gPC2+gPC3+gPC4+gPC5+gPC6+ePC1+ePC2+ePC3+ePC4+ePC5,
                                  family = "poisson", nAGQ=0, data= mod_data, control = glmerControl(optimizer = "nloptwrap"))
    	    model_lrt <- anova(null_model, full_model)
            res = data.frame("SNP" = snp, "GENE" = gene,
                 "BETA" = summary(full_model)$coefficients[2,][1], #G beta
                "SE" = summary(full_model)$coefficients[2,][2], #G se
                "P" = model_lrt$`Pr(>Chisq)`[2])
            all_res = rbind(all_res,res)}, error=function(cond){return(NA)})
    }
    write.csv(all_res, paste0(args$res_filepath, args$celltype_geno,'_',args$lead_snp,'_csaQTL_test_',args$celltype_expr,'_sceQTLs_',gene,'.csv'), row.names = FALSE)
}

