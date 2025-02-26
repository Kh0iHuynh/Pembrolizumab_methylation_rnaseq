# Load necessary libraries
library(dplyr)

# Step 1: Read the first file (file1.csv) with gene expression data
file1 <-read.table("deseq_result_cohortb_survival.txt",header = TRUE,sep = "\t",row.names =1)

file2 <-read.table("mann_whitney_methylation_cohortb_result.txt",header = TRUE,sep = " ")


print("done load files")



library(biomaRt)


# Convert rownames (Ensemble gene names) into a column 'gene'
file1$gene <- rownames(file1)

# Set up connection to Ensembl via biomaRt
# Query Ensembl's gene annotation data for chromosome, start, and end positions
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Extract chromosome, start, and end positions for each gene
gene_info <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'),
                   filters = 'ensembl_gene_id',
                   values = file1$gene, # Use Ensembl gene IDs from file1
                   mart = ensembl)


print("done extract info")
# Merge the extracted data (gene_info) with file1 based on gene IDs
file1 <- merge(file1, gene_info, by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)

# Add 5000 base pairs to the start and end positions
file1 <- file1 %>%
  mutate(
    start_hg38 = start_position - 5000,   # Subtract 5000 from start
       chromosome_name = paste0("chr", chromosome_name),  # Add "chr" prefix to chromosome_name
    end_hg38 = end_position + 5000       # Add 5000 to end
  )


# Filter rows from file2 based on overlap with regions from file1
# Create an empty data frame to store matched rows from file2
matched_rows <- data.frame()


print("start matching file")

# Loop through each row of file2 and compare with file1 regions
for (i in 1:nrow(file2)) {
  # Extract chr, start, and end for the current row of file2
  chr2 <- file2$chr_hg38[i]
  start2 <- file2$start_hg38[i]
  end2 <- file2$end_hg38[i]
  
  # Filter file1 to find rows that match the chromosome and overlap with the region in file2
  overlaps <- file1 %>%
    filter(chromosome_name == chr2 &  # Match chromosome
           start_hg38 <= end2 &  # Start of file1 region is before or at the end of file2 region
           end_hg38 >= start2)   # End of file1 region is after or at the start of file2 region
  
  # If there is any overlap, add the row from file2 to matched_rows
  if (nrow(overlaps) > 0) {
    matched_rows <- rbind(matched_rows, file2[i, , drop = FALSE])
  }
}

# Write the matched rows to a new file (matched_file.csv)
write.table(matched_rows, file = "matched_file.txt", row.names = FALSE, sep = "\t", quote = FALSE)
# Print a message confirming the result
cat("Matched rows have been written to 'matched_file.csv'.\n")

