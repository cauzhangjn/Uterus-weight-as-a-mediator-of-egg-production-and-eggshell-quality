#!/bin/bash
## Author: Zhang junnan
## SBATCH -c 20
## Step 2：Calculate PEER factors and pcs, peer in another server analysis
name=$1
module load EIGENSOFT/7.2.1
result=xx/xx/${name}
bin=xx/xx/code
pcs_res=xx/xx/${name}/${name}_pcs
cd ${result}
mkdir ${name}_pcs
python ${bin}/compute_genotype_pcs.py ${result}/${name}.vcf.gz --keep -o ${pcs_res}
