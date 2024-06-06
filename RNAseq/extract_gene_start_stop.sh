#!/bin/bash
#SBATCH --account=modrek_1129
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --array=1


while read line; do
        grep -e "$line" ../reference_genome/gencode.v45_gene_annotation_table.txt | awk -v gene=$line '{print gene, $3,$4,$5}' >> DE_gene_withstartstop.txt
done < DE_gene.txt
