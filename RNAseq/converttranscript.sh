#!/bin/bash

#####
# convert trasncript name to only gene name
#####
cut -f1 -d" " deseq_result_cohort.txt | cut -f1 -d"." > delete.txt
cut -f2,3,4,5,6,7 -d" " deseq_result_cohort.txt > delete2.txt
paste -d" " delete.txt delete2.txt > deseq_result_cohort_ensgene.txt



cut -f1 -d" " deseq_result_survival.txt | cut -f1 -d"." > delete.txt
cut -f2,3,4,5,6,7 -d" " deseq_result_survival.txt > delete2.txt
paste -d" " delete.txt delete2.txt > deseq_result_survival_ensgene.txt


cut -f1 -d" " deseq_result_interaction.txt | cut -f1 -d"." > delete.txt
cut -f2,3,4,5,6,7 -d" " deseq_result_interaction.txt > delete2.txt
paste -d" " delete.txt delete2.txt > deseq_result_interaction_ensgene.txt
