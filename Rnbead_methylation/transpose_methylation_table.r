library(tidyverse)

#####
# turn of scientific notation
#####
options(scipen = 999)

####
# read in normalized count
# convert SRR col name to row values
####
samplecount <- as.data.frame(read.csv('hg38.methylation_beta_value.txt', header = TRUE, sep = '\t'))


samplecount_long <- samplecount %>% pivot_longer(cols = DF2CAB_tumor_gDNA_2:UT34SJ_tumor_gDNA_, values_to = "methylation_beta", names_to = "Sample_ID")

samplecount_long[1:4,1:3]
#####
# read in metadata
#####


# join metadata to methylation beta data
#####

samples <- read.csv('deseq_metadata_withouthyphen.csv', header = TRUE, sep = ' ')

#####
# join metadata to normalized count data
#####

alldata = samplecount_long %>% full_join(samples,join_by(Sample_ID))

write.table(alldata, file = "transposed_methylation.txt", col.names = TRUE , row.names = FALSE, quote =FALSE,sep ="\t")
