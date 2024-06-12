#!/bin/bash
#SBATCH --account=modrek_1129
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --array=1

#####
# intialize shell to run conda
# activate nanopore conda environment
#####
eval "$(conda shell.bash hook)"
conda activate scseq


module load samtools/1.17 bedtools2/2.30.0
module load  gcc/11.3.0  openblas/0.3.20
module load r/4.3.2 bowtie2/2.4.2

module load pkgconf libpng libjpeg-turbo libtiff zlib bzip2 libxml2 curl fontconfig freetype graphite2 pcre harfbuzz fribidi glib openssl icu4c automake autoconf libtool m4 xz libiconv

module load trimmomatic/0.39 fastqc/0.11.9
module load bedtools2/2.30.0


Rscript methylresolver.r
#Rscript methylation_anova.r
