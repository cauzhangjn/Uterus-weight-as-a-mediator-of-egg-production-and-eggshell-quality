#!/bin/bash
#BATCH -p Cnode
#SBATCH -c 10
#SBATCH --mem=10G

bedtools intersect -a UW100_filtered.bed -b chicken.bed -wa -wb > UW100_genes_intersect.bed
