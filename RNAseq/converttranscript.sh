#!/bin/bash

#####
# convert trasncript name to only gene name
#####
cut -f1 -d" " anova_normalized_count_result.txt | cut -f1 -d"." > delete.txt
cut -f2,3,4,5 -d" " anova_normalized_count_result.txt > delete2.txt
paste -d" " delete.txt delete2.txt > deseq_anova_result.txt
