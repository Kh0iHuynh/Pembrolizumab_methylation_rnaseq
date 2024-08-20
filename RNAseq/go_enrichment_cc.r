library(org.Hs.eg.db)
library(GO.db)
library(GOstats)
library(tidyverse)


#####
# turn of scientific notation
#####
options(scipen = 999)

#####
# survival go_analysis
#####

survival_data1 <- read.table("deseq/deseq_significant_survival.txt",header = TRUE, sep=" ")

entrez_id <- read.table("ensembletoentrez.txt",header = TRUE, sep=" ")

survival_data =survival_data1 %>% left_join(entrez_id,join_by(ensgene))

###convert DE gene ensg to entrez
#data <- read.table("DE_gene.txt",header = TRUE)
#entrez_id = mapIds(org.Hs.eg.db, keys = data$ensgene, keytype="ENSEMBL", column = "ENTREZID")
#entrez_id[1]
#write.table(entrez_id,file = "delete.txt", quote = FALSE)
#transposed_entrez_id <-as.data.frame(t(entrez_id))



geneup = unique(survival_data[survival_data$log2FoldChange > 1,"entrez_id"])
genedown = unique(survival_data[survival_data$log2FoldChange < -1,"entrez_id"])
universegene = unique(survival_data$entrez_id)

goupparam = new("GOHyperGParams", 
	geneIds = geneup,
	universeGeneIds = universegene , 
	annotation = "org.Hs.eg.db",
	ontology = "CC",
	pvalueCutoff= 0.05, 
	conditional=FALSE, 
	testDirection="over")


godownparam = new("GOHyperGParams",
        geneIds = genedown,
        universeGeneIds = universegene ,
        annotation = "org.Hs.eg.db",
        ontology = "CC",
        pvalueCutoff= 0.05,
        conditional=FALSE,
        testDirection="over")

up_go_run =hyperGTest(goupparam)
down_go_run =hyperGTest(godownparam)
summary_up = summary(up_go_run)
summary_down= summary(down_go_run)
write.table(summary_up,file = "goanalysis_cc_survival_up.txt",row.names=FALSE, quote = FALSE,sep="\t")
write.table(summary_down,file = "goanalysis_cc_survival_down.txt", row.names=FALSE, quote = FALSE,sep="\t")


#####
# response go_analysis
#####


response_data1 <- read.table("deseq/deseq_significant_cohort.txt",header = TRUE, sep=" ")

entrez_id <- read.table("ensembletoentrez.txt",header = TRUE, sep=" ")

response_data =response_data1 %>% left_join(entrez_id,join_by(ensgene))
#entrez_id = mapIds(org.Hs.eg.db, keys = response_data$ensgene, keytype="ENSEMBL", column = "ENTREZID")
#entrez_id[1]
#write.table(entrez_id,file = "delete.txt", quote = FALSE)
#transposed_entrez_id <-as.data.frame(t(entrez_id))


annotation= loadDb("gencode.v32.annotation.TxDb")

geneup = unique(response_data[response_data$log2FoldChange > 1,"entrez_id"])
genedown = unique(response_data[response_data$log2FoldChange < -1,"entrez_id"])
universegene = unique(response_data$entrez_id)

goupparam = new("GOHyperGParams", 
	geneIds = geneup,
	universeGeneIds = universegene , 
	annotation = "org.Hs.eg.db",
	ontology = "CC",
	pvalueCutoff= 0.05, 
	conditional=FALSE, 
	testDirection="over")


godownparam = new("GOHyperGParams",
        geneIds = genedown,
        universeGeneIds = universegene ,
        annotation = "org.Hs.eg.db",
        ontology = "CC",
        pvalueCutoff= 0.05,
        conditional=FALSE,
        testDirection="over")

up_go_run =hyperGTest(goupparam)
down_go_run =hyperGTest(godownparam)
summary_up = summary(up_go_run)
summary_down= summary(down_go_run)
write.table(summary_up,file = "goanalysis_cc_cohort_up.txt",row.names=FALSE, quote = FALSE,sep="\t")
write.table(summary_down,file = "goanalysis_cc_cohort_down.txt", row.names=FALSE, quote = FALSE,sep="\t")

#####
# interaction go_analysis
#####

interaction_data1 <- read.table("deseq/deseq_significant_interaction.txt",header = TRUE, sep=" ")

entrez_id <- read.table("ensembletoentrez.txt",header = TRUE, sep=" ")

interaction_data =interaction_data1 %>% left_join(entrez_id,join_by(ensgene))
#entrez_id = mapIds(org.Hs.eg.db, keys = interaction_data$ensgene, keytype="ENSEMBL", column = "ENTREZID")
#entrez_id[1]
#write.table(entrez_id,file = "delete.txt", quote = FALSE)
#transposed_entrez_id <-as.data.frame(t(entrez_id))


annotation= loadDb("gencode.v32.annotation.TxDb")

geneup = unique(interaction_data[interaction_data$log2FoldChange > 1,"entrez_id"])
genedown = unique(interaction_data[interaction_data$log2FoldChange < -1,"entrez_id"])
universegene = unique(interaction_data$entrez_id)

goupparam = new("GOHyperGParams", 
	geneIds = geneup,
	universeGeneIds = universegene , 
	annotation = "org.Hs.eg.db",
	ontology = "CC",
	pvalueCutoff= 0.05, 
	conditional=FALSE, 
	testDirection="over")


godownparam = new("GOHyperGParams",
        geneIds = genedown,
        universeGeneIds = universegene ,
        annotation = "org.Hs.eg.db",
        ontology = "CC",
        pvalueCutoff= 0.05,
        conditional=FALSE,
        testDirection="over")

up_go_run =hyperGTest(goupparam)
down_go_run =hyperGTest(godownparam)
summary_up = summary(up_go_run)
summary_down= summary(down_go_run)
write.table(summary_up,file = "goanalysis_cc_interaction_up.txt",row.names=FALSE, quote = FALSE,sep="\t")
write.table(summary_down,file = "goanalysis_cc_interaction_down.txt", row.names=FALSE, quote = FALSE,sep="\t")



