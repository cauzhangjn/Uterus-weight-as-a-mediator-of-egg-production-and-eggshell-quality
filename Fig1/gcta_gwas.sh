#!/bin/bash
#BATCH -p Cnode
#SBATCH -c 10
#SBATCH --mem=10G


trait=$1  
type=$2 
threads=$3

cd workspace
# mkdir -p $trait
cd $trait
#### 0. get duli_snp_number
plink --bfile ../$type --chr-set 39 --chr 1-39 --allow-extra-chr --indep-pairwise 50 10 0.2 --make-founders
wc -l plink.prune.in > duli_snp_number.txt

#### 1. get pheno
j=`awk -v col=$trait 'NR==1 {for (i=1; i<=NF; i++) if ($i == col) print i}' ../${type}.txt`
awk '{print $1,$2,$'$j'}' ../${type}.txt > ${trait}.txt

#### 2. build geno matrix
## build G-matrix
gcta64 --bfile ../${type} --make-grm  --autosome-num 39 --out ${trait} --thread-num ${threads}
gcta64 --grm ${trait} --make-bK-sparse 0.05 --out ${trait}_sp_grm

## calculate pca
gcta64 --grm ${trait} --pca 10 --thread-num ${threads} --out ${trait}
awk '{NF=7}1' ${trait}.eigenvec > ${trait}.pc.txt
Rscript  pca.R ${trait}

#### 3. calculate h2
gcta64 --grm ${trait} --pheno ${trait}.txt --reml --out ${trait} --thread-num ${threads}

#### 4. run gwas
# --covar fixed.txt
gcta64 --fastGWA-mlm --bfile ../${type} --grm-sparse ${trait}_sp_grm --pheno ${trait}.txt --qcovar ${trait}.pc.txt --out ${trait} --thread-num ${threads} --autosome-num 39 

#### 5. gwas plot 
Rscript  gwas_plot.R ${trait}
#mv Rect_Manhtn.P-value.jpg  $trait.Manhattan.jpg
#mv QQplot.P-value.jpg  $trait.QQplot.jpg
