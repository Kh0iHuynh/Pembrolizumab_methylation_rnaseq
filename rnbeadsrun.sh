#!/bin/bash

#SBATCH --account=modrek_1129
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00          

module load  gcc/11.3.0  openblas/0.3.20
module load r/4.3.2
module load conda/23.10.0
module load zip/3.0
#R also need impute package

#intialize shell to run conda
eval "$(conda shell.bash hook)"
conda activate rnbeads

Rscript rnbead.r

