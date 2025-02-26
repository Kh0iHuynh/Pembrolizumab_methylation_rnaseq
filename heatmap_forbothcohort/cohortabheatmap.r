# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)

pdf("combined_cohortab_heatmap.pdf", width = 7, height = 7)
# Load the data (replace 'file1.csv' and 'file2.csv' with your actual file paths)
file11 <- read.table("normalized_genebothset.txt", header = TRUE, sep = " ")  # Adjust path and separator if needed
#file 1 header :ensgene SRR_ID normalized_count Sample_ID Sentrix_ID Sentrix_Position Sex Age Overall_Survival RANO Relapses_Count OS_Status IFNG_response KPS Group Cohort


file4 <- read.table("genebothset.txt", header = TRUE, sep = " ")  # Adjust path and separator if needed

# file 2 heade : baseMean        log2FoldChange  lfcSE   stat    pvalue  padj


# Remove any numbers and the period following 'ensgene' using regex
file11$ensgene <- gsub("\\.\\d+$", "", file11$ensgene)
file4$ensgene <- gsub("\\.\\d+$", "", file4$ensgene)
head(file11)
head(file4)

samplecondition = c("A","B")
genecondition = c("A","B")


for (i in 1:length(samplecondition)) {
file1 = file11 %>% filter(Cohort == samplecondition[i])

for (j in 1:length(genecondition)) {

# Sorting file2 by p_value to get the top IDs with the lowest p_value
top_ids_file2 <- file4 %>% filter(Cohort == genecondition[j])  # Keep top 100 IDs with lowest p-value
# Filtering file1 based on the selected IDs from file2
filtered_file1 <- file1[file1$ensgene %in% top_ids_file2$ensgene, ]
##################
# Reshape the data to have Sample_ID as column names and ID as row names
reshaped_data <- filtered_file1 %>%
  select(ensgene, Sample_ID, normalized_count) %>%  # Select relevant columns
  spread(key = Sample_ID, value = normalized_count)  # Reshape the data

# Set ID as rownames
rownames(reshaped_data) <- reshaped_data$ensgene
reshaped_data <- reshaped_data[, -1]  # Remove the ID column as it is now the row names


################33
# extract factors
#################

file3 <- read.table("combined_deseq_metadata.csv", header = TRUE, sep = " ")

file2 = file3 %>% filter(Cohort == samplecondition[i])


# Assuming reshaped_data and file2 are already loaded
# Extract column names from reshaped_data (these are the Sample_IDs)
sample_ids_reshaped <- colnames(reshaped_data)

# Extract the Sample_ID column from file2 (adjust column names if needed)
sample_ids_file2 <- file2$Sample_ID

# Step 1: Match Sample_IDs in reshaped_data with those in file2
# Ensure that only Sample_IDs in reshaped_data are retained in file2
matched_data <- file2 %>% 
  filter(Sample_ID %in% sample_ids_reshaped) %>%
  arrange(match(Sample_ID, sample_ids_reshaped))  # Ensure the order matches that in reshaped_data


##########33
# make sure reshaped data has same columns
###########3

# Step 1: Extract Sample_IDs from both reshaped_data and matched_data
sample_ids_reshaped <- colnames(reshaped_data)
sample_ids_matched <- matched_data$Sample_ID

# Step 2: Filter columns in reshaped_data to only keep those in sample_ids_matched
reshaped_data_filtered <- reshaped_data[, sample_ids_reshaped %in% sample_ids_matched]

# Step 3: Reorder the columns of reshaped_data to match the order of Sample_IDs in matched_data
reshaped_data_filtered <- reshaped_data_filtered[, match(sample_ids_matched, colnames(reshaped_data_filtered))]

# Step 2: Extract Overall_Survival and Cohort from matched file2
overall_survival_list <- matched_data$Overall_Survival
cohort_list <- matched_data$Cohort


# Optional: Create named lists for further use
list_overall_survival <- overall_survival_list
list_cohort <- cohort_list

##############33
# zcale and plot
###############
print("start plotting")

#zscale row wise

scaled_data <- scale(t(reshaped_data_filtered))

###replace value to keep constant -1 1 range
scaled_data[scaled_data > 1 ] <- 1
scaled_data[scaled_data < -1] <- -1

#mat_z <- t(scale(t(reshaped_data_filtered)))

mat_z <- t(scaled_data)


#col annotation
mat_z[1,]
rownames(mat_z)


Overall_Survival = overall_survival_list

dim(data.frame(Overall_Survival))
print(mat_z)


ha = HeatmapAnnotation(Overall_Survival = Overall_Survival,col = list(Overall_Survival = c("Over" =  "red", "Under" = "blue")))


plot <- Heatmap(mat_z, name = "Z-score", show_row_names = TRUE, show_column_names = TRUE,top_annotation = ha,column_split = data.frame(Overall_Survival),column_title = paste0("Cohort ",samplecondition[i]," Gene Expression for \n Cohort ",genecondition[j]," DE gene by Survival"),row_names_gp = gpar(fontsize = 8, rot = 45),cluster_columns=FALSE)
print(plot)
}
}




#####################
# CPG
######################
file11 <- read.table("../transposed_methylation.txt", header = TRUE, sep = "\t")  # Adjust path and separator if needed

file4 <- read.table("cpgbothset.txt", header = TRUE, sep = " ")  # Adjust path and separator if needed


samplecondition = c("A","B")
genecondition = c("A","B")


for (i in 1:length(samplecondition)) {
file1 = file11 %>% filter(Cohort == samplecondition[i])

for (j in 1:length(genecondition)) {
# Sorting file2 by p_value to get the top IDs with the lowest p_value

top_ids_file2 = file4 %>% filter(Cohort == genecondition[j])


# Filtering file1 based on the selected IDs from file2
filtered_file1 <- file1[file1$ID %in% top_ids_file2$ID, ]
print("filtered_file1")
head(filtered_file1)
##################
# Reshape the data to have Sample_ID as column names and ID as row names
reshaped_data <- filtered_file1 %>%
  select(ID, Sample_ID, methylation_beta) %>%  # Select relevant columns
  spread(key = Sample_ID, value = methylation_beta)  # Reshape the data

# Set ID as rownames
rownames(reshaped_data) <- reshaped_data$ID
reshaped_data <- reshaped_data[, -1]  # Remove the ID column as it is now the row names

print("reshaped")
head(reshaped_data)



################33
# extract factors
#################

file2 <- read.table("combined_deseq_metadata_withouthyphen.csv", header = TRUE, sep = " ")
# Assuming reshaped_data and file2 are already loaded
# Extract column names from reshaped_data (these are the Sample_IDs)
sample_ids_reshaped <- colnames(reshaped_data)

# Extract the Sample_ID column from file2 (adjust column names if needed)
sample_ids_file2 <- file2$Sample_ID

# Step 1: Match Sample_IDs in reshaped_data with those in file2
# Ensure that only Sample_IDs in reshaped_data are retained in file2
matched_data <- file2 %>%
  filter(Sample_ID %in% sample_ids_reshaped) %>%
  arrange(match(Sample_ID, sample_ids_reshaped))  # Ensure the order matches that in reshaped_data


##########33
# make sure reshaped data has same columns
###########3

# Step 1: Extract Sample_IDs from both reshaped_data and matched_data
sample_ids_reshaped <- colnames(reshaped_data)
sample_ids_matched <- matched_data$Sample_ID

# Step 2: Filter columns in reshaped_data to only keep those in sample_ids_matched
reshaped_data_filtered <- reshaped_data[, sample_ids_reshaped %in% sample_ids_matched]

# Step 3: Reorder the columns of reshaped_data to match the order of Sample_IDs in matched_data
reshaped_data_filtered <- reshaped_data_filtered[, match(sample_ids_matched, colnames(reshaped_data_filtered))]




# Step 2: Extract Overall_Survival and Cohort from matched file2
overall_survival_list <- matched_data$Overall_Survival
cohort_list <- matched_data$Cohort

# Step 3: Print the lists to check if they are aligned correctly
print("Overall Survival List:")
print(overall_survival_list)

print("Cohort List:")
print(cohort_list)

# Optional: Create named lists for further use
list_overall_survival <- overall_survival_list
list_cohort <- cohort_list

##############33
# zcale and plot
###############
print("start plotting")

#zscale row wise

scaled_data <- scale(t(reshaped_data_filtered))

###replace value to keep constant -1 1 range
scaled_data[scaled_data > 1 ] <- 1
scaled_data[scaled_data < -1] <- -1

#mat_z <- t(scale(t(reshaped_data_filtered)))

mat_z <- t(scaled_data)


#col annotation
mat_z[1,]
rownames(mat_z)


Overall_Survival = overall_survival_list

dim(data.frame(Overall_Survival))
dim(mat_z)


ha = HeatmapAnnotation(Overall_Survival = Overall_Survival,col = list(Overall_Survival = c("Over" =  "red", "Under" = "blue")))


plot <- Heatmap(mat_z, name = "Z-score", show_row_names = TRUE, show_column_names = TRUE,top_annotation = ha,column_split = data.frame(Overall_Survival),column_title =  paste0("Cohort ",samplecondition[i]," CpG methylation for \n Cohort ",genecondition[j]," DM CpG by Survival"),row_names_gp = gpar(fontsize = 8, rot = 45),cluster_columns=FALSE)
print(plot)


}
}


dev.off()
