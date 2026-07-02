#!/bin/bash
#SBATCH -p thunder
#SBATCH -c 15
#SBATCH --mem 15gb
#SBATCH -J twas


python ${SPrediXcan} \
--gwas_file ${gwas_file} \
--snp_column SNP  --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column pval_cor  \
--model_db_path ${db}  \
--covariance ${covariances} \
--keep_non_rsid  --model_db_snp_key rsid \
--throw \
--output_file $output/spredixcan/chicken_${tissue}_PM_spredixcan.csv

###############

awk -F',' '$5 < 0.05' $output/spredixcan/chicken_${tissue}_PM_spredixcan.csv > $output/spredixcan/filtered_chicken_${tissue}_PM_spredixcan.csv
