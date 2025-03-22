library(scales)
library(fgsea)
library(msigdbr)
library(ReactomePA)
library(enrichplot)
library(cowplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(EnhancedVolcano)
library(AnnotationHub)
library(multtest)
library(SeuratDisk)
library(patchwork)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(tidyverse)


pdf("goanalysis/goanalysisgse_degene.pdf", width = 7, height = 7)

args = commandArgs(trailingOnly=TRUE)

options(future.globals.maxSize = 8000 * 1024^2)
options(scipen = 999)


# Configure the Ensembl mart connection for human genes

########
# create gene annotation
#########
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
#####
# load in seurat object from h5Seurat
####
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


#=------------
# convert to hgnc symbol
#------------
markers1 <- read.table("deseq_result.txt",header = TRUE, sep = "\t")
markers = markers1 %>% dplyr::filter(pvalue < 0.05)
# we want the log2 fold change
# Convert markers to a data frame if not already
df <- as.data.frame(markers)
head(df)
# Save avg_log2FC as the original gene list for reference
original_gene_list <- df$log2FoldChange
# Retrieve mapping of hgnc_symbol to entrezgene_id using biomaRt
mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = df$ensgene,
  mart = mart)
head(mapping)
head(df)
df$ensembl_gene_id <- df$ensgene
# Merge mapping into the original data frame to include entrezgene_id
df <- merge(df, mapping, by = "ensembl_gene_id", all.x = TRUE)
# Remove rows with missing entrezgene_id
df <- na.omit(df)
df <- df[!is.na(df$hgnc_symbol) & df$hgnc_symbol != "", ]
head(df)
# Create a named vector for gene list using avg_log2FC as values and entrezgene_id as names
gene_list <- df$log2FoldChange
names(gene_list) <- df$hgnc_symbol
# Remove NA values if any remain (as a safety measure)
gene_list <- na.omit(gene_list)
# Sort the list in decreasing order (required by clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

head(gene_list)

#########33
# process loop
############



library(msigdbr)
library(dplyr)

# Corrected msigdb_categories list
msigdb_categories <- list(
  H = NULL,  # H has no subcategories
  C2 = c("CGP","CP:BIOCARTA","CP:PID","CP:REACTOME", "CP:WIKIPATHWAYS"),
  C3 = c("MIR:MIRDB","TFT:GTRD"),
  C4 = c("CGN","CM"),
  C5 = c("GO:BP","GO:CC", "GO:MF","HPO"),
  C6 = NULL,  # C6 has no subcategories
  C7 = c("IMMUNESIGDB","VAX"),
  C8 = NULL  # C8 has no subcategories
)





print("going into loop")

# Loop through the main categories and subcategories
for (main_category in names(msigdb_categories)) {
  cat("Main Category:", main_category, "\n")

  # Check if the category has subcategories
  if (!is.null(msigdb_categories[[main_category]])) {
    # Loop through the subcategories for each main category
    for (subcategory in msigdb_categories[[main_category]]) {
      cat("  Subcategory:", subcategory, "\n")

      # Dynamically process based on the category and subcategory variables
      # Retrieve gene sets based on main category and subcategory dynamically
      category_sets <- msigdbr(species = "Homo sapiens", category = main_category, subcategory = subcategory)
      category_sets_list <- category_sets %>% 
                            split(.$gs_name) %>% 
                            lapply(function(x) x$gene_symbol)
      cat("  Gene sets retrieved for", main_category, ":", subcategory, "\n")

      # Prepare ranked gene list based on the aggregated values
      #gene_list <- markers$log2FoldChange
      #names(gene_list) <- rownames(markers)

      # Check if gene_list is empty
      if (length(gene_list) == 0) {
        cat("Warning: gene_list is empty for cluster:", cluster, "and condition:", conditions[[j]], "\n")
        next  # Skip this iteration if gene_list is empty
      }

      # Check overlap between gene list and pathways
      pathway_genes <- unique(unlist(category_sets_list))
      gene_overlap <- intersect(names(gene_list), pathway_genes)
      cat("Number of overlapping genes:", length(gene_overlap), "\n")

      # Perform GSEA
      fgsea_results <- fgsea(pathways = category_sets_list,
                             stats = gene_list,
                             minSize = 0,  # Minimum size of pathways to consider
                             maxSize = 8000)  # Maximum size of pathways to consider

      # Check if fgsea_results is empty before printing
      if (nrow(fgsea_results) > 0) {
        fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ", "))
        print(head(fgsea_results))
        write.table(fgsea_results, file = paste0("goanalysis/gsea_", main_category, "_", subcategory, "_for_degene.txt"),
                    sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

        ###########
        # start plotting
        ###########
        # Assuming fgsea_results is the data frame with your results
        # Filter for pathways with p-value less than 0.05
        fgsea_significant <- fgsea_results %>%
          dplyr::filter(pval < 0.05)

        # Check if there are any significant pathways
        if (nrow(fgsea_significant) > 0) {
          # Select top 3 and bottom 3 by NES
          top_bottom_3 <- fgsea_significant %>%
            arrange(desc(NES)) %>%  # Sort by NES descending
            slice_head(n = 3) %>%
            bind_rows(fgsea_significant %>%
                        arrange(NES) %>%  # Sort by NES ascending
                        slice_head(n = 3))

          # You can adjust the size to be based on the number of genes in the leading edge or another metric
          top_bottom_3 <- top_bottom_3 %>%
            mutate(size = sapply(leadingEdge, function(x) length(strsplit(x, ",")[[1]])))  # Size based on the number of genes in leadingEdge

          # Now create the lollipop plot
          plot <- ggplot(top_bottom_3, aes(x = size, y = reorder(pathway, NES), fill = NES)) +
            geom_bar(stat = "identity", width = 0.7) +  # Bar plot with size on the x-axis
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.5, 2.5)) +  # Blue to white to red gradient with midpoint at 0
            ggtitle(paste0(main_category," " ,subcategory, " for de gene")) +
            theme_classic() +
            labs(x = "Size",
                 y = "Pathway",
                 fill = "NES") +
            theme(legend.position = "right",
                  axis.text.y = element_text(hjust = 0),  # Align pathway names to the right
                  axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate size labels for clarity
            scale_y_discrete(labels = function(x) str_wrap(gsub("_", " ", x), width = 20))  # Replace _ with space and wrap text

          ggsave(filename =  paste0("goanalysis/gsea_", main_category, "_", subcategory, "_for_degene.pdf"), plot = plot, width = 7, height = 7, dpi = 600)

          print(plot)
        } else {
          cat("No significant pathways (p-value < 0.05) found.\n")
        }

        ###########
        # Done plotting
        #############
      }
    }
  } else {
    # If no subcategories, retrieve the gene sets directly for the category
    cat("  No subcategories\n")

    # Get gene sets for the main category (e.g., C6, C8)
    main_category_sets <- msigdbr(species = "Homo sapiens", category = main_category)
    main_category_sets_list <- main_category_sets %>%
                               split(.$gs_name) %>%
                               lapply(function(x) x$gene_symbol)
    cat("  Gene sets retrieved for category:", main_category, "\n")

    # Prepare ranked gene list based on the aggregated values
    #gene_list <- markers$log2FoldChange
    #names(gene_list) <- rownames(markers)

    # Check if gene_list is empty
    if (length(gene_list) == 0) {
      cat("Warning: gene_list is empty for cluster:", cluster, "and condition:", conditions[[j]], "\n")
      next  # Skip this iteration if gene_list is empty
    }

    # Check overlap between gene list and pathways
    pathway_genes <- unique(unlist(main_category_sets_list))
    gene_overlap <- intersect(names(gene_list), pathway_genes)
    cat("Number of overlapping genes:", length(gene_overlap), "\n")

    # Perform GSEA
    fgsea_results <- fgsea(pathways = main_category_sets_list,
                           stats = gene_list,
                           minSize = 0,  # Minimum size of pathways to consider
                           maxSize = 8000)  # Maximum size of pathways to consider

    # Check if fgsea_results is empty before printing
    if (nrow(fgsea_results) > 0) {
      fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ", "))
      print(head(fgsea_results))
      write.table(fgsea_results, file = paste0("goanalysis/scseq2_gsea_", main_category, "_for_degene.txt"),
                  sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

      ###########
      # start plotting
      ###########
      # Assuming fgsea_results is the data frame with your results
      # Filter for pathways with p-value less than 0.05
      fgsea_significant <- fgsea_results %>%
        dplyr::filter(pval < 0.05)

      # Check if there are any significant pathways
      if (nrow(fgsea_significant) > 0) {
        # Select top 3 and bottom 3 by NES
        top_bottom_3 <- fgsea_significant %>%
          arrange(desc(NES)) %>%  # Sort by NES descending
          slice_head(n = 3) %>%
          bind_rows(fgsea_significant %>%
                      arrange(NES) %>%  # Sort by NES ascending
                      slice_head(n = 3))

        # You can adjust the size to be based on the number of genes in the leading edge or another metric
        top_bottom_3 <- top_bottom_3 %>%
          mutate(size = sapply(leadingEdge, function(x) length(strsplit(x, ",")[[1]])))  # Size based on the number of genes in leadingEdge

        # Now create the lollipop plot
        plot <- ggplot(top_bottom_3, aes(x = size, y = reorder(pathway, NES), fill = NES)) +
          geom_bar(stat = "identity", width = 0.7) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.5, 2.5)) +
          ggtitle(paste0(main_category, " for de gene")) +
          theme_classic() +
          labs(x = "Size", y = "Pathway", fill = "NES") +
          theme(legend.position = "right",
                axis.text.y = element_text(hjust = 0),
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_y_discrete(labels = function(x) str_wrap(gsub("_", " ", x), width = 20))  # Replace _ with space and wrap text

        ggsave(filename = paste0("goanalysis/gsea_", main_category, "_for_degene.pdf"), plot = plot, width = 7, height = 7, dpi = 600)
        print(plot)
      } else {
        cat("No significant pathways (p-value < 0.05) found.\n")
      }

      ###########
      # Done plotting
      #############
    }
  }

  cat("\n")
}
dev.off()
