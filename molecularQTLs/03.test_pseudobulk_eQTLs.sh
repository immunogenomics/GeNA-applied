#!/bin/bash
debug=1 #un-commenting this line will cause script to execute commands rather than submitting                                        

# Test for associations of csaQTL lead SNPs to pseudobulked expression within each major cell type

csaQTL_celltypes=( "Myeloid" "NK" "NK" "NK" "NK")
lead_snps=( "15:80263217:C:T" "2:111851212:C:T" "11:128070535:A:G" "12:10583611:C:T" "19:16441973:G:A" )
chrs=( 15 2 11 12 19 )
n_snps=${#lead_snps[@]}

for i in $(eval echo "{0..$(($n_snps-1))}")
do
    for expr_celltype in "allcells"
    do
	celltype=${csaQTL_celltypes[i]}
	lead_snp=${lead_snps[i]}
	chr=${chrs[i]}
	
	out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/pseudobulk/${celltype}_${lead_snp}_csaQTL_test_${expr_celltype}_eQTLs"
	geno_doses="/data/srlab/lrumker/MCSC_Project/cna-qtl/results/geno_munge/cis_snps/${celltype}_${lead_snp}_cis.DS.vcf.gz"
	geno_sample_ids="/data/srlab/lrumker/datasets/onek1k/geno/sample_labels/sample_list_chr${chr}.txt"
	peer_resid_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/pseudobulk/inputs/${expr_celltype}_PEER_residuals.csv"
	genes_file="/data/srlab/lrumker/MCSC_Project/cna-qtl/mask_cis/results/cis_genes_${celltype}_${lead_snp}.tsv"
	
	command="Rscript test_pseudobulk_eQTLs.R --genes_file $genes_file --out_path $out_path \
                     --peer_resid_file $peer_resid_file --geno_sample_ids $geno_sample_ids --geno_doses $geno_doses"
	echo $command
	
	if [ -z "$debug" ]
	then
            mkdir -p out
            name="pseudobulk_norm_${celltype}"
            bsub -J $name -q big \
		-oo out/${name}.out \
		-eo out/${name}.err \
		-R 'select[hname!=cn001]' \
		-R 'select[hname!=cn002]' \
		-R 'select[hname!=cn003]' \
		-R 'select[hname!=cn004]' \
		-R 'select[hname!=cn005]' \
		-R 'select[hname!=cn007]' \
		"$command"
	else
            LSB_JOBINDEX=1
            eval $command
            #exit
	fi
    done
done