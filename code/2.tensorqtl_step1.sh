#!/bin/bash
## Author: Zhang junnan
## SBATCH -c 5
## Step 1 : Generate normalized expression 
## Officially start analysis
## Step 1 : Generate normalized expression in BED format ##Generate ${type}.expression.bed.gz

for type in chuiti huichang luanchao pengda shierzhichang xiaqiunao zigong ganzhang mangchang yixian
do
    path=xx/xx
    result=xx/xx/${type}
    bin=xx/xx/code
    data=xx/xx/${type}
    cd ${result}
    python ${bin}/eqtl_prepare_expression.py ${data}/${type}_TPM.txt ${data}/${type}_count.txt ${path}/chicken.gtf \
        ${result}/${type}_sample_id.txt ${result}/${type}_chr_list ${result}/${type} \
        --tpm_threshold 0.1 \
        --count_threshold 6 \
        --sample_frac_threshold 0.2 \
        --normalization_method tmm
done
