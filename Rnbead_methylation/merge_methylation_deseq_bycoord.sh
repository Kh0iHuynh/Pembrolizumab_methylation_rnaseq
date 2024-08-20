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


####
# convert space to tab for DE_gene_withstartstop.txt
# convert gene start and stop to 10k region 
# start - 5000 , stop + 5000
####

awk '{print $1,$2,$3 - 5000,$4 + 5000}' deseq/DE_gene_withstartstop.txt > DE_gene_temp.txt

#####
# 
#####

awk '{print $2,$3,$4,$1,$5,$6,$7}' methylation/anova_methylation_result.txt | tail -n +2 | sort -k1,1 -k2,2n -k3,3n | sed 's/ /\t/g' > methylationresulttemp.txt



#####
# extract start stop and gene name
#####

echo "chromosome        start   stop    ensgene baseMean        log2FoldChange  lfcSE   stat    pvalue  padj    DE_factor	chr_cpg	start_cpg	end_cpg	cpg_ID	betalog10survival	betalog10cohort	betalog10survivalbycohort" > de_and_methylation_anova_result.txt


awk 'NR==FNR{a[$1]=$0;next}{if($1 in a){print a[$1],$0}}' DE_gene_temp.txt deseq/deseq_significant_survival.txt | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,"survival"}' | sed 's/ /\t/g' > deseq_significant_temp.txt


awk 'NR==FNR{a[$1]=$0;next}{if($1 in a){print a[$1],$0}}' DE_gene_temp.txt deseq/deseq_significant_cohort.txt | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,"cohort"}' | sed 's/ /\t/g' >> deseq_significant_temp.txt


awk 'NR==FNR{a[$1]=$0;next}{if($1 in a){print a[$1],$0}}' DE_gene_temp.txt deseq/deseq_significant_interaction.txt | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,"interaction"}' | sed 's/ /\t/g' >> deseq_significant_temp.txt


grep -v -e "chrM" deseq_significant_temp.txt | sort -k1,1 -k2,2n -k3,3n > delete.txt
mv delete.txt deseq_significant_temp.txt
#####
# combine methylation anova and deseq result
#####

bedtools intersect -wa -wb \
    -a deseq_significant_temp.txt \
    -b methylationresulttemp.txt >> de_and_methylation_anova_result.txt
