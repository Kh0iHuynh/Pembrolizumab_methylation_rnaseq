# Load required libraries
library(tidyverse)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(variancePartition)
library(BiocParallel)
library(ggplot2)
library(pheatmap)

# Turn off scientific notation
options(scipen = 999)

# 1. Read in Data and Metadata
samples <- read.csv('deseq_metadata.csv', header = TRUE, sep = ' ')
rownames(samples) <- samples$SRR_ID
samples

files <- file.path("/project/modrek_1129/KhoiHuynh/chow_rna", "salmon", samples$SRR_ID, "editted.quant.sf")
names(files) <- samples$SRR_ID

# 2. Create tx2gene Table
txdb <- makeTxDbFromGFF(file = "/project/modrek_1129/KhoiHuynh/reference_genome/gencode.v45.annotation.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# 3. Import Quant Files with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi$counts[is.na(txi$counts)] <- 0

# 4. Construct DESeq2 Object and Normalize Data
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Overall_Survival)
dds <- DESeq(ddsTxi)

# Results for Overall Survival
res <- results(dds, contrast = c("Overall_Survival", "Over", "Under"))
res <- res[order(res$padj),]

# Clean rownames of DESeq2 results
rownames(res) <- gsub("\\..*", "", rownames(res))  # Remove '.' and anything after
write.table(res, file = "deseq_result_cohortb_survival.txt", col.names = TRUE, quote = FALSE, row.names = TRUE,sep = "\t")

# Plot MA-plot
plotMA(res, ylim = c(-5, 5), main = "Survival")

# Obtain normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Clean rownames of normalized counts
rownames(normalized_counts) <- gsub("\\..*", "", rownames(normalized_counts))  # Remove '.' and anything after
write.table(normalized_counts, file = "normalized_counts_cohortb.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


# 5. Variance Partitioning with Sex and Age
# Align sample names
colnames(normalized_counts) <- gsub("\\..*", "", colnames(normalized_counts)) # Remove extensions if any
samples <- samples[rownames(samples) %in% colnames(normalized_counts), ]
normalized_counts <- normalized_counts[, rownames(samples)]

# Verify alignment
stopifnot(all(rownames(samples) == colnames(normalized_counts)))

# Prepare ExpressionSet object for variancePartition
expr_matrix <- as.matrix(normalized_counts)
pheno <- samples
expr_obj <- ExpressionSet(assayData = expr_matrix, phenoData = AnnotatedDataFrame(pheno))

# Define the formula for variance partitioning
form <- ~ Overall_Survival + Age + Sex
# Use Multicore for Speed (Optional)
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# Run variancePartition
varPart <- fitExtractVarPartModel(expr_obj, form, pheno, BPPARAM = param)

# Sort and plot variance partition results
vp <- sortCols(varPart)
plotVarPart(vp)

# Clean rownames of variance partitioning results
rownames(vp) <- gsub("\\..*", "", rownames(vp))  # Remove '.' and anything after
write.table(vp, file = "percent_variance_explained_cohortb.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 6. Summary Output
cat("DESeq2 analysis completed.\n")
cat("Variance Partitioning results saved to 'percent_variance_explained.txt'.\n")

data <- read.table("percent_variance_explained_cohortb.txt",header =  TRUE, row.names = 1,sep ="\t")
data2 <- read.table("deseq_result_cohortb_survival.txt",header =  TRUE, row.names = 1,sep ="\t")

data2 = data2 %>% mutate(gene=rownames(data2))
data = data %>% mutate(gene=rownames(data))

data3 <- inner_join(data, data2, by = c("gene" = "gene"))
write.table(data3, file = "joint_deseq_result_percent_variance_explained_cohortb.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Completion Message
cat("Joint DESeq2 and variance partition results saved.\n")
