library( "DESeq2" )
library(ggplot2)
library(GenomicFeatures)
library("pheatmap")
#####
# read in data and metadata file
# data file = "genename SRR1_tpm SRR2_tpm"
# meta data = "ID (SRRID), factor1, factor2"
#####
#countData <- as.matrix(read.csv('deseq_data.txt', header = TRUE,sep = ' ',row.names=1))
samples <- read.csv('deseq_metadata.csv', header = TRUE, sep = ' ')
rownames(samples) <- samples$SRR_ID
samples$SRR_ID

files <- file.path("/project/modrek_1129/KhoiHuynh/chow_rna","salmon", samples$SRR_ID, "editted.quant.sf")
names(files) <- samples$SRR_ID
files
#####
# create tx2gene
#####
txdb <- makeTxDbFromGFF(file="/project/modrek_1129/KhoiHuynh/reference_genome/gencode.v45.annotation.gtf.gz")
saveDb(x=txdb, file = "gencode.v45.annotation.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME") 
head(tx2gene)
#####
# import salmon quant files
#####
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

head(txi$counts)
#####
# set any reference transcript to 0 if NA
#####
txi$counts[is.na(txi$counts)] <- 0

#####
# turn of scientific notation
#####
options(scipen = 999)
#####
# construct deseq object
# model : overall survival
#####

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Overall_Survival)

dds <- DESeq(ddsTxi)

#####
# assign result to res variable
# 
#####i

res <- results( dds, contrast = c("Overall_Survival", "Over", "Under" ) )
res <- res[order(res$padj),]
write.table(res, file = "deseq_result_survival.txt", col.names= TRUE, quote = FALSE)
# plot MA-plot
plotMA(res, ylim = c(-5, 5),main= " Survival")

#####
# obtaining normalized count
#####

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="deseq_survival_normalized_counts.txt", sep="\t", quote=F, col.names=TRUE)


#####
# Cohort model 
#####
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Cohort)

dds <- DESeq(ddsTxi)

#####
# assign result to res variable
# sort result by pvalue
#####
res <- results( dds, contrast = c("Cohort", "A", "B" ) )
res <- res[order(res$padj),]
write.table(res, file = "deseq_result_cohort.txt", col.names= TRUE, quote = FALSE)
plotMA(res, ylim = c(-5, 5),main= "Cohort")

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="deseq_cohort_normalized_counts.txt", sep="\t", quote=F, col.names=TRUE)

#####
# IFN-G_response and Overall_survival
#####
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Overall_Survival + Cohort + Overall_Survival:Cohort)

dds <- DESeq(ddsTxi)

#####
# assign result to res variable
# sort result by pvalue
#####
res <- results(dds)
res <- res[order(res$padj),]
write.table(res, file = "deseq_result_interaction.txt", col.names= TRUE, quote = FALSE)
plotMA(res, ylim = c(-5, 5),main= "Interaction")

