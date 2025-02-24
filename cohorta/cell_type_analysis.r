library(dplyr)
library(dplyr)
library(tidyverse)
#####
# turn of scientific notation
#####
options(scipen = 999)

#######
# tranform methylresolver result
#######
resolverdata <- as.data.frame(read.csv('MethylResolver.txt', header = TRUE,row.names = 1, sep = '\t'))
resolverdata$Sample_ID <- row.names(resolverdata)

#resolverdata2 = resolverdata %>% select(Sample_ID,abs_Mon, abs_Dendritic, abs_Macro, abs_Neu, abs_Eos, abs_Treg, abs_Tnaive, abs_Tmem, abs_CD8, abs_NK, abs_Bcell, Purity)

resolverdata2 = resolverdata %>% select(Sample_ID,Mon, Dendritic,Macro, Neu, Eos, Treg, Tnaive, Tmem, CD8, NK, Bcell)

samplecount_long <- resolverdata2 %>% pivot_longer(cols = Mon:Bcell , values_to = "Frequency", names_to = "Cell_type")

samplecount_long[1:3,]

write.table(samplecount_long, file = "methylation_celltype_count.txt", ,append= FALSE,quote = FALSE,row.names=FALSE)
#####
#######


samplemetadata <- read.table("deseq_metadata_withouthyphen.csv",header = TRUE, sep=" ")

type_count <- read.table("methylation_celltype_count.txt",header = TRUE,sep =" ")

type_count3 = type_count %>% mutate(Percentage = Frequency*100)

immune_content2 = type_count3 %>% left_join(y=samplemetadata,join_by(Sample_ID)) %>% group_by(Sample_ID) %>% summarise( total_immune = sum(Percentage[Cell_type != "Purity"]))  

immune_content3 = type_count3 %>% left_join(y=samplemetadata,join_by(Sample_ID))
immune_content  = immune_content2 %>% left_join(y=samplemetadata,join_by(Sample_ID))

immune_content2
immune_content
na.omit(immune_content)
####
# anova by cell_type
####

ANOVA = immune_content3 %>% group_by(Cell_type) %>%
    do(out = anova(lm(Percentage~Overall_Survival,na.action=na.omit,data=.))) %>%
    mutate(-pf(out[1,3]/out[2,3],out[1,1],out[2,1],lower.tail = FALSE, log.p = TRUE)/log(10)) %>%
    select(-out)

na.omit(ANOVA)
ANOVA

write.table(ANOVA,file = "anova_by_cell_type_result.txt",append= FALSE,quote = FALSE,row.names=FALSE)

####
# anova by total immune cell count 
####

anova(lm(total_immune~Overall_Survival,data=immune_content))

