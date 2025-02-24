library(tidyverse)

#####
# turn off scientific notation
#####
options(scipen = 999)

####
# Read in normalized count
# Convert SRR column name to row values
####

# Read in the methylation data
samplecount <- as.data.frame(read.csv('../transposed_methylation.txt', header = TRUE, sep = '\t'))

# Prepare data by log-transformation (adding 5 to avoid log(0))
alldata2 <- samplecount %>% mutate(transform_beta = log(methylation_beta + 5))

# Perform Mann-Whitney U test (Wilcoxon test) for each combination of ID, chr, and region
mann_whitney_results <- alldata2 %>%
  filter(Cohort == "B") %>%  # Filter for Cohort B
  group_by(ID, chr_hg38, start_hg38, end_hg38) %>%
  do({
    # Perform the Mann-Whitney U test on Overall_Survival and transform_beta
    test_result <- wilcox.test(transform_beta ~ Overall_Survival, data = ., na.action = na.omit)
    # Return p-value
    data.frame(p_value = test_result$p.value)
  })

# Extract and view the results
mann_whitney_results <- mann_whitney_results %>%
  mutate(log10_pvalue = -log10(p_value))  # Convert p-values to -log10(p-value)

# Write the result to a file
write.table(mann_whitney_results, file = "mann_whitney_methylation_cohortb_result.txt", append = FALSE, quote = FALSE, row.names = FALSE)

# View the results (optional)
print(mann_whitney_results)

