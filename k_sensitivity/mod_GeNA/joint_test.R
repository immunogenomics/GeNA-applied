library(argparse)
library(Rmpfr)
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--chisq_per_nampc_file",type="character")
parser$add_argument("--ks_file",type="character")
parser$add_argument("--outfile",type="character")
args <- parser$parse_args()
#cat("\n\n****")
#print(args)
#cat("****\n\n")

all_res = read.table(args$chisq_per_nampc_file, header=FALSE)
all_res = all_res**2 # T-squared
ks = read.table(args$ks_file, header=FALSE)[,1]

mht_correct <-function(raw_p, len_ks, prec_bits = 100){asNumeric(mpfr(1, prec_bits)-(mpfr(1, prec_bits)-mpfr(raw_p, prec_bits))**len_ks)}

for(k in ks){
    all_res[paste0("k",k,"_P")] = apply(as.matrix(rowSums(all_res[,c(1:k), drop=FALSE]), ncol=1), 1, pchisq, df = k, lower = F)
}

# In this modified script, P values for all tested k are exported
write.table(all_res[,paste0("k",ks,"_P"), drop=FALSE], paste0(args$outfile, "_all_k"), quote=FALSE, row.names=FALSE, sep = "\t")
all_res['P'] = apply(as.matrix(all_res[,paste0("k",ks,"_P"), drop=FALSE], ncol=length(ks)), 1, min)

# In this modified script, we are evaluating the result for all values of k *individually* so 
# we don't correct across the number of k values considered
#all_res['P'] = apply(as.matrix(all_res[,"P", drop=FALSE]), 1, mht_correct, length(ks))

all_res['k']= apply(as.matrix(all_res[,paste0("k",ks,"_P"), drop=FALSE], ncol=length(ks)), 1, which.min)
all_res['k'] = ks[all_res$k]
write.table(all_res[,c("P", "k")], args$outfile, quote=FALSE, row.names=FALSE, sep = "\t")