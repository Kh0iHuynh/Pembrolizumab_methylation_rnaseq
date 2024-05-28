#!/bin/bash
#SBATCH --account=modrek_1129
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --time=48:00:00
#SBATCH --array=1-23

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
#####
# Set environment variable
#####

array_namelist="./namelist.txt"
input_dir="./sra"
input_trim="./trimmed"
output_dir="./salmon"
qc_dir="./qc_dir"

input_namelist=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $array_namelist)
echo "$input_namelist"

#mkdir -p $input_trim
#mkdir -p $qc_dir/raw
#mkdir -p $qc_dir/trim
#mkdir -p $output_dir/$input_namelist

#####
# fastqc raw reads
#####

#mkdir -p $qc_dir/raw/$input_namelist
#fastqc $input_dir/${input_namelist}.fastq.gz -o $qc_dir/raw/$input_namelist

#####
# trimm read 
#Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
#Remove leading low quality or N bases (below quality 3) (LEADING:3)
#Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
#Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
#Drop reads below the 36 bases long (MINLEN:36)
#####

#trimmomatic SE -threads 8 -phred33 $input_dir/${input_namelist}.fastq.gz $input_trim/${input_namelist}.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#####
# fastqc trimmed reads
#####

#mkdir -p $qc_dir/trim/$input_namelist
#fastqc $input_trim/${input_namelist}.fastq.gz -o $qc_dir/trim/$input_namelist

#####
# alignment using salmon
#####

#salmon quant -i /project/modrek_1129/KhoiHuynh/reference_genome/genecode_v45_hg38 -l A -r $input_trim/${input_namelist}.fastq.gz -p 8 --validateMappings -o $output_dir/$input_namelist


#####
# make quant file for deseq
#####
echo "Name	Length	EffectiveLength	TPM	NumReads" > $output_dir/${input_namelist}/editted.quant.sf
awk -v sample=$input_namelist '{print $0,sample}' $output_dir/$input_namelist/quant.sf | sed 's/|/ /g' | tail -n +2 | awk '{print $1,$9,$10,$11,$12}' | sort -k1,1 | sed 's/ /\t/g' > $output_dir/${input_namelist}/editted.quant.sf.1

# compare and output NA to transcript not in quant
awk 'NR==FNR{a[$1]=$0;next}{if($1 in a){print a[$1]}else{print $1,"NA","NA","NA","NA"}}' $output_dir/${input_namelist}/editted.quant.sf.1 /project/modrek_1129/KhoiHuynh/reference_genome/gencode.v45_transcript.txt | grep -v "Geneid" | sed 's/ /\t/g' > $output_dir/${input_namelist}/editted.quant.sf.2

cat $output_dir/${input_namelist}/editted.quant.sf.1 $output_dir/${input_namelist}/editted.quant.sf.2 >> $output_dir/${input_namelist}/editted.quant.sf

#####
# make tpm quant with ENSG instead
#####
echo "Name	Length	EffectiveLength	TPM	NumReads" > $output_dir/${input_namelist}/ENSG.editted.quant.sf

awk -v sample=$input_namelist '{print $0,sample}' $output_dir/$input_namelist/quant.sf | sed 's/|/ /g' | tail -n +2 | awk '{print $2,$9,$10,$11,$12}' | sort -k1,1 | sed 's/ /\t/g' >> $output_dir/${input_namelist}/ENSG.editted.quant.sf

#####
# make tpm file with all unique ENSG to prepare for TPM matrix
# for normalization
#####
echo "$input_namelist" > $output_dir/${input_namelist}.uniquetranscript.txt
awk 'NR==FNR{a[$1]=$4;next}{if($1 in a){print $1,a[$1]}else{print $1,"0"}}' $output_dir/${input_namelist}/ENSG.editted.quant.sf uniquegene.txt | grep -v "Geneid" | sed 's/ /\t/g' | sort -k1,1 | awk '{print $2}' >> $output_dir/${input_namelist}.uniquetranscript.txt

