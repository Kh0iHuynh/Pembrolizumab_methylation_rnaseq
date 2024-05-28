library(tidyverse)

#####
# turn of scientific notation
#####
options(scipen = 999)

####
# read in normalized count
# convert SRR col name to row values
####
samplecount <- as.data.frame(read.csv('normalized_counts.txt', header = TRUE, sep = '\t',row.names =1))

samplecount$ensgene <-row.names(samplecount)

samplecount_long <- samplecount %>% pivot_longer(cols = SRR8112646:SRR8112674, values_to = "normalized_count", names_to = "SRR_ID")

samplecount_long[1:4,1:3]
#####
# read in metadata
#####


samples <- read.csv('deseq_metadata.csv', header = TRUE, sep = ' ')

#####
# join metadata to normalized count data
#####

alldata = samplecount_long %>% full_join(samples,join_by(SRR_ID))
alldata[1:3,]

write.table(alldata, file = "normalized_count_with_metadata.txt", col.names = TRUE , row.names = FALSE, quote =FALSE)

#####
# ANOVA
# alldata header:
# ensgene SRR_ID normalized_count Sample_ID Sentrix_ID Sentrix_Position Sex Age Overall_Survival RANO Relapses_Count OS_Status IFNG_response KPS
#####

alldata2 = alldata %>% mutate(transform_count = log(normalized_count +5 ))
ANOVA = alldata2 %>% group_by(ensgene) %>% 
    do(out = anova(lm(transform_count~Overall_Survival + IFNG_response + Overall_Survival:IFNG_response,na.action=na.omit,data=.))) %>%
    mutate(countlog10survival = -pf(out[1,3]/out[3,3],out[1,1],out[3,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    mutate(countlog10response = -pf(out[2,3]/out[3,3],out[2,1],out[3,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    mutate(countlog10survivalbyresponse = -pf(out[3,3]/out[4,3],out[3,1],out[4,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    mutate(residualcount = out[4,3]) %>%
    select(-out)

na.omit(ANOVA)

write.table(ANOVA,file = "anova_normalized_count_result.txt",append= FALSE,quote = FALSE,row.names=FALSE)
