#!/bin/bash
#SBATCH -p Cnode
#SBATCH -c 10
#SBATCH --mem=50G
tissue=$1

workpath=xx/smr
outpath=xx/smr/output
genotype_plink=xx/all
smr_type_file=/smr_prepare_file

cd ${workpath}


./smr_Linux --bfile ${genotype_plink} \
--gwas-summary UW100.ma \
--beqtl-summary ${smr_type_file}/eQTL_${tissue}_merged.cis_qtl_pairs \
--heidi-mtd 1 \
--out ${outpath}/eQTL_${tissue}
