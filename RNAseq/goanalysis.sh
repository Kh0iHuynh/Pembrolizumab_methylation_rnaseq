#!/bin/bash
#SBATCH --account=modrek_1129
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00


module load samtools/1.17 bedtools2/2.30.0
module load  gcc/11.3.0  openblas/0.3.20
module load r/4.3.2

module load pkgconf libpng libjpeg-turbo libtiff zlib bzip2 libxml2 curl fontconfig freetype graphite2 pcre harfbuzz fribidi glib openssl icu4c automake autoconf libtool m4 xz libiconv
#####


head -n 1 deseq_significant_survival.txt > survival_up.txt
head -n 1 deseq_significant_survival.txt > survival_down.txt
head -n 1 deseq_significant_survival.txt > survival_neutral.txt

head -n 1 deseq_significant_cohort.txt > cohort_up.txt
head -n 1 deseq_significant_cohort.txt > cohort_down.txt
head -n 1 deseq_significant_cohort.txt > cohort_neutral.txt


head -n 1 deseq_significant_interaction.txt > interaction_up.txt
head -n 1 deseq_significant_interaction.txt > interaction_down.txt
head -n 1 deseq_significant_interaction.txt > interaction_neutral.txt


awk '{if( $3 >= -1 && $3 <= 1 ) print $0}' deseq_significant_survival.txt | grep -v -e "ensgene" >> survival_neutral.txt
awk '{if(  $3 > 1 ) print $0}' deseq_significant_survival.txt | grep -v -e "ensgene" >> survival_up.txt

awk '{if(  $3 < -1 ) print $0}' deseq_significant_survival.txt | grep -v -e "ensgene" >> survival_down.txt

awk '{if( $3 >= -1 && $3 <= 1 ) print $0}' deseq_significant_cohort.txt | grep -v -e "ensgene" >> cohort_neutral.txt
awk '{if(  $3 > 1 ) print $0}' deseq_significant_cohort.txt | grep -v -e "ensgene" >> cohort_up.txt

awk '{if(  $3 < -1 ) print $0}' deseq_significant_cohort.txt | grep -v -e "ensgene" >> cohort_down.txt

awk '{if( $3 >= -1 && $3 <= 1 ) print $0}' deseq_significant_interaction.txt | grep -v -e "ensgene" >> interaction_neutral.txt
awk '{if(  $3 > 1 ) print $0}' deseq_significant_interaction.txt | grep -v -e "ensgene" >> interaction_up.txt

awk '{if(  $3 < -1 ) print $0}' deseq_significant_interaction.txt | grep -v -e "ensgene" >> interaction_down.txt
