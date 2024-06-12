library(tidyverse)
library(lme4)
#####
# turn of scientific notation
#####
options(scipen = 999)

####
# read in normalized count with metadata
# file header:
# chr_hg38        start_hg38      end_hg38        chr_hg19        start_hg19      end_hg19        ID      Sample_ID       methylation_beta        SRR_ID  Sentrix_ID      Sentrix_Position        Sex     Age     Overall_Survival        RANO    Relapses_Count  OS_Status       Cohort   KPS     Group
####
samplecount2 <- as.data.frame(read.csv('transposed_methylation.txt', header = TRUE, sep = '\t'))

samplecount = samplecount2 %>% mutate(transform_count = log(methylation_beta +5 ))

lme4sex = samplecount %>% group_by(ID) %>% do(Model = as.data.frame(VarCorr(lmer(transform_count ~ (1|Sex) + (1 | Sex:Overall_Survival) + (1|Overall_Survival) + (1|Sex:Cohort) + (1|Overall_Survival:Cohort:Sex) + (1|Group) + (1|Sex:Group), data = .)))) %>% unnest(Model)


#####
# variance explained by Sex 
# for Overall_Survival
#####
PvarSexandSurvival <- function(x)
{
100*(x$vcov[x$grp=="Sex"]/(x$vcov[x$grp=="Sex"]+x$vcov[x$grp=="Sex:Overall_Survival"]))
}

PvarS = lme4sex %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarSexandSurvival(.))))

write.table(PvarS,file="Pvar.sexandsurvival.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("sex by survival done")
#####
# variance explained by Sex
# for INFG_Cohort
#####
PvarSexandResponse <- function(x)
{
100*(x$vcov[x$grp=="Sex"]/(x$vcov[x$grp=="Sex"]+x$vcov[x$grp=="Sex:Cohort"]))
}

PvarS = lme4sex %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarSexandResponse(.))))

write.table(PvarS,file="Pvar.sexandCohort.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("Sex by Cohort done")


######
# Sex by group
######

PvarSexandGroup <- function(x)
{
100*(x$vcov[x$grp=="Sex"]/(x$vcov[x$grp=="Sex"]+x$vcov[x$grp=="Sex:Group"]))
}

PvarS = lme4sex %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarSexandGroup(.))))

write.table(PvarS,file="Pvar.sexandgroup.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("Sex by group done")


#####
# Variance explained by Sex 
# for Survival and Response interaction
#####

PvarSinteraction <- function(x)
{
100*(x$vcov[x$grp=="Sex:Cohort"]/(x$vcov[x$grp=="Sex:Overall_Survival"]+x$vcov[x$grp=="Overall_Survival:Cohort:Sex"]+ x$vcov[x$grp=="Sex:Cohort"]))
}

PvarS = lme4sex %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarSinteraction(.))))

write.table(PvarS,file="Pvar.sexandinteraction.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("Sex by Interaction done")
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################


lme4Age = samplecount %>% group_by(ID) %>% do(Model = as.data.frame(VarCorr(lmer(transform_count ~ (1|Age) + (1 | Age:Overall_Survival) + (1|Overall_Survival) + (1|Age:Cohort) + (1|Overall_Survival:Cohort:Age) + (1|Group) + (1|Age:Group), data = .)))) %>% unnest(Model)


#####
# variance explained by Age
# for Overall_Survival
#####
PvarAgeandSurvival <- function(x)
{
100*(x$vcov[x$grp=="Age"]/(x$vcov[x$grp=="Age"]+x$vcov[x$grp=="Age:Overall_Survival"]))
}

PvarS = lme4Age %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarAgeandSurvival(.))))

write.table(PvarS,file="Pvar.ageandsurvival.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("Age by survival done")
#####
# variance explained by Age
# for INFG_Cohort
#####
PvarAgeandResponse <- function(x)
{
100*(x$vcov[x$grp=="Age"]/(x$vcov[x$grp=="Age"]+x$vcov[x$grp=="Age:Cohort"]))
}

PvarS = lme4Age %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarAgeandResponse(.))))

write.table(PvarS,file="Pvar.ageandCohort.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("Age by Cohort done")
#####
# Variance explained by Age
# for Survival and Response interaction
#####

PvarAgeinteraction <- function(x)
{
100*(x$vcov[x$grp=="Age:Cohort"]/(x$vcov[x$grp=="Age:Overall_Survival"]+x$vcov[x$grp=="Overall_Survival:Cohort:Age"] +x$vcov[x$grp=="Age:Cohort"]))
}

PvarS = lme4Age %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarAgeinteraction(.))))

write.table(PvarS,file="Pvar.ageandinteraction.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("age by interaction done")


#####
# variance explained by Age
# for Group
#####
PvarAgeandGroup <- function(x)
{
100*(x$vcov[x$grp=="Age"]/(x$vcov[x$grp=="Age"]+x$vcov[x$grp=="Age:Group"]))
}

PvarS = lme4Age %>% group_by(ID) %>% do(data.frame(PvarS = round(PvarAgeandGroup(.))))

write.table(PvarS,file="Pvar.ageandgroup.txt",row.names=FALSE,sep="\t",quote=FALSE)

print("Age by group done")



########################################################
########################################################





