#!/bin/bash

# Generates PEER factors and residuals for gene expression

debug=1 #un-commenting this line will cause script to execute commands rather than submitting                                        

# conda activate peer
# PEER installation troubleshooting help: https://github.com/mz2/peer/issues/9

for celltype in "Myeloid" "NK"
do
    out_path="/data/srlab/lrumker/MCSC_Project/cna-qtl/eqtls/results/trans_eQTL_candidates/pseudobulk/${celltype}"
    expr_file="${out_path}_samplesXgenes_norm_invnt.csv"
    meta_file="${out_path}_samples_meta.csv"
    covs_list="age,sex,gPC1,gPC2,gPC3,gPC4,gPC5,gPC6"
    n_peer_factors=20
    
    command="Rscript run_PEER.R --expr_file $expr_file --n_peer_factors $n_peer_factors \
              --meta_file $meta_file  --covs_list $covs_list --out_path $out_path"
    echo $command
    
    if [ -z "$debug" ]
    then
        mkdir -p out
        name="PEER_${celltype}"
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
            exit
    fi
done