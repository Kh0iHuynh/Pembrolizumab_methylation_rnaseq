#!/bin/bash
#SBATCH --account=modrek_1129
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --array=1


module load  gcc/11.3.0  openblas/0.3.20
module load r/4.3.2

#Rscript delete.r
#Rscript methylation_anova.r
Rscript lme4.r
