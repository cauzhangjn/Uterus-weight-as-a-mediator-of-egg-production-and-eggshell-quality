#!/bin/bash
## Author: Zhang junnan
## SBATCH -c 20
## Step 3：Combine PEER factors and pcs
bin=/xx/xx/code
cd ${bin}
Rscript ${bin}/combine_peer_and_pc.R
