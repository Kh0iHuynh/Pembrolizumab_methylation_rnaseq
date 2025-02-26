pdf("venndiagram.pdf", width = 7, height = 7)

data <-read.table("bothab.txt",header = TRUE,sep = " ")

# bothab.txt is a text file contain gene_name and Cohort for all DE genes from cohort a and b
# Create a list of unique genes by cohort
gene_sets <- split(data$gene_name, data$Cohort)
gene_sets <- lapply(gene_sets, unique)
# Load the VennDiagram package
library(VennDiagram)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = names(gene_sets),  # Names of the cohorts
  filename = NULL,  # Don't save to a file just yet
  output = TRUE
)

# Display the Venn diagram
grid.draw(venn.plot)

# Save the Venn diagram as a PNG file
venn.diagram(
  x = gene_sets,
  category.names = names(gene_sets),
  filename = "venn_diagram.png",
  output = TRUE
)

dev.off()
