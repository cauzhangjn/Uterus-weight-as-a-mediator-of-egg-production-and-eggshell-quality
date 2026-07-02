#!/bin/bash
## Author: Zhang junnan


## Step 4 : Run Tensorqtl (Website: https://github.com/broadinstitute/tensorqtl)
prefix=$1
bin=/xx/xx/code
genoetype=/xx/xx/${prefix}/${prefix}
phenotype=/xx/xx/${prefix}/${prefix}.expression.bed.gz
covariates=/xx/xx/${prefix}/${prefix}.combined_covariates.txt
cd /xx/xx/${prefix}
mkdir ${prefix}_qtl
outpath=/xx/xx/${prefix}/${prefix}_qtl

#cd /xx/xx/${prefix}

# 1) cis-QTL mapping: compute cis nominal associations for all variant-gene pairs
python -m tensorqtl \
           ${prefix} ${phenotype} ${prefix} \
           --mode cis_nominal \
           --covariates ${covariates} \
           --maf_threshold 0.05 \
           -o ${outpath}

# 2) cis-QTL mapping: compute cis permutations and define eGenes
python -m tensorqtl \
           ${prefix} ${phenotype} ${prefix} \
           --mode cis \
           --covariates ${covariates} \
           --maf_threshold 0.05 \
           -o ${outpath}

# 3) perform fdr test
cd ${outpath}
Rscript ${bin}/tensorqtl_fdr.R ${outpath}/${prefix}.cis_qtl.txt.gz ${outpath}/${prefix}.cis_qtl.txt 0.05
cat ${outpath}/${prefix}.cis_qtl.txt | awk '{if (NR==1 || $20 =="TRUE") print}'  > ${outpath}/${prefix}.cis_qtl.txt.fdr.txt

cd /xx/xx/${prefix}
# 4) cis-QTL mapping: conditionally independent QTLs
cis_output=${outpath}/${prefix}.cis_qtl.txt.fdr.txt
python -m tensorqtl \
            ${prefix} ${phenotype} ${prefix} \
            --mode cis_independent \
            --covariates ${covariates} \
            --cis_output ${cis_output} \
            --maf_threshold 0.05
# 5) trans-QTL mapping
python -m tensorqtl \
           ${prefix} ${phenotype} ${prefix} \
           --covariates ${covariates} \
           --mode trans \
           --maf_threshold 0.05 --output_text 

echo "Done!"
