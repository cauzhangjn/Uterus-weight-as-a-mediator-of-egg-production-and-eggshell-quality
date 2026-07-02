#!/bin/bash
## Author: Zhang junnan
## SBATCH -c 5
## Step 2：Calculate PEER factors and pcs
name=$1
bin=xx/code
outpath=xx/${name}
cd ${outpath}
module load R/4.0.2

## Run

Rscript ${bin}/PEER.R ${outpath}/${name}.expression.bed.gz ${name} 30 --output_dir ${outpath}
