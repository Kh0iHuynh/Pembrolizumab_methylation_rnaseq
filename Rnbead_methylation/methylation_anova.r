###########################################
# SURVIVAL
###########################################

library(tidyverse)

#####
# turn of scientific notation
#####
options(scipen = 999)

####
# read in normalized count
# convert SRR col name to row values
####

# chr_hg38        start_hg38      end_hg38        chr_hg19        start_hg19      end_hg19 ID      Strand  Sample_ID       methylation_beta        SRR_ID  Sentrix_ID       Sentrix_Position        Sex     Age     Overall_Survival        RANO    Relapses_Count   OS_Status       IFNG_response   KPS     Group   Cohort

samplecount <- as.data.frame(read.csv('transposed_methylation.txt', header = TRUE, sep = '\t'))


alldata2 = samplecount %>% mutate(transform_beta = log(methylation_beta +5 ))
ANOVA = alldata2 %>% group_by(ID,chr_hg38,start_hg38,end_hg38) %>% 
    do(out = anova(lm(transform_beta~Overall_Survival + Cohort + Overall_Survival:Cohort,na.action=na.omit,data=.))) %>%
    mutate(betalog10survival = -pf(out[1,3]/out[3,3],out[1,1],out[3,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    mutate(betalog10cohort = -pf(out[2,3]/out[3,3],out[2,1],out[3,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    mutate(betalog10survivalbycohort = -pf(out[3,3]/out[4,3],out[3,1],out[4,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    mutate(residualcount = out[4,3]) %>%
    select(-out)

na.omit(ANOVA)

write.table(ANOVA,file = "anova_methylation_result.txt",append= FALSE,quote = FALSE,row.names=FALSE)
###########################################
# RESPONSE
###########################################
