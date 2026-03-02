# Load packages
library(BiocManager)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(DOSE)
library(tidyverse)
library("DESeq2")
library("tximport")
library(ggplot2)
library("apeglm")
library("ashr")
library(pheatmap)
#1. Data Acquisition----
#Load metadata
sampleTable <- read.csv("stage_Metadata.csv",row.names=1)
#Assigning file path and names
files <- file.path("Gene count tables",paste0(rownames(sampleTable),".fastq.genes.results"))
names(files) <- rownames(sampleTable)
#Import 
#Gene counts already mapped to reference with rsem tx2gene will not be used 
txi <- tximport(files, type="rsem")
#Sanity checks
dim(txi$counts)
head(txi$counts)
#Checking that txi data and metadata have same names
all.equal(colnames(txi$counts), rownames(sampleTable))

#2. Deseq DGE Analysis----
#Early biofilm as reference
levels(factor(sampleTable$Stage))
#Deseq analysis
dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~Stage)
dds <- DESeq(dds)

resultsNames(dds)
res_mature_vs_early <- results(dds, name="Stage_Mature.Biofilm_vs_Early.biofilm")
res_thin_vs_early   <- results(dds, name="Stage_Thin.Biofilm_vs_Early.biofilm")
res_mature_vs_thin <- results(dds, contrast=c("Stage", "Mature Biofilm", "Thin Biofilm"))
summary(res_mature_vs_early)
summary(res_mature_vs_thin)
summary(res_thin_vs_early)

#Likelihood ratio test to see how much stage improve the model
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
#63% of genes significantly affected by stage 
summary(res_lrt)

#3. Deseq plots----
# Shrinkage and plots
#Shrinkage using lfcshrink
mature_vs_earlyLFC <- lfcShrink(dds, coef="Stage_Mature.Biofilm_vs_Early.biofilm", type="apeglm")

Thin_vs_earlyLFC <- lfcShrink(dds, coef="Stage_Thin.Biofilm_vs_Early.biofilm", type="apeglm")
mature_vs_ThinLFC <- lfcShrink(dds, contrast=c("Stage", "Mature Biofilm", "Thin Biofilm"), type="ashr")

#MA plots
plotMA(mature_vs_earlyLFC, ylim=c(-2,2))
plotMA(Thin_vs_earlyLFC, ylim=c(-2,2))
plotMA(mature_vs_ThinLFC, ylim=c(-2,2))

#Create data frames for each pairwise results
Mature_vs_Early_df <- as.data.frame(mature_vs_earlyLFC)
Thin_vs_Early_df <- as.data.frame(Thin_vs_earlyLFC)
Mature_vs_Thin_df <- as.data.frame(mature_vs_ThinLFC)
makeVolcano <- function(res_df) {
  
  # Capture dataframe name
  df_name <- deparse(substitute(res_df))
  clean_name <- gsub("_", " ", sub("_df$", "", df_name))
  
  # Add gene column
  res_df$gene <- rownames(res_df)
  
  # Remove NA values first
  res_df <- na.omit(res_df)
  
  # Define significance
  res_df$significant <- ifelse(
    res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
    ifelse(res_df$log2FoldChange > 0, "Up", "Down"),
    "Not Sig"
  )
  
  # Plot
  ggplot(res_df,
         aes(x = log2FoldChange,
             y = -log10(pvalue),
             color = significant)) +
    geom_point() +
    scale_color_manual(values = c("Down" = "blue",
                                  "Not Sig" = "gray",
                                  "Up" = "red")) +
    labs(x = "Log2 Fold Change",
         y = "-Log10 p-value",
         title = paste("Volcano Plot:", clean_name)) +
    theme_minimal() +
    theme(legend.position = "right")
}
makeVolcano(Mature_vs_Early_df)
makeVolcano(Mature_vs_Thin_df)
makeVolcano(Thin_vs_Early_df)

#Heatmap
#Using likelihood ratio test results to represent all stages as a time course not just pairwise
# Remove NA values
res_lrt <- na.omit(res_lrt)

# Select top 20 genes by p values from non-NA genes
top_genes <- head(order(res_lrt$padj), 20)
gene_names <- rownames(res_lrt)[top_genes]

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)
# Store counts in a matrix for the heatmap
mat <- assay(vsd)[gene_names, ]

# Add annotation for Sample ID  and Stage
annotation_df <- sampleTable[, c("Sample.ID", "Stage")]
colnames(annotation_df) <- c("Sample ID", "Stage")

# Create heatmap
pheatmap(mat, 
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = FALSE,
)

#PCA Plot

# Same VST used above for the heatmap
vsd <- vst(dds)

# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = c("Stage"), returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# GGplot code to display Stage by colour
ggplot(pca_data, aes(x = PC1, y = PC2, color = Stage)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") + coord_fixed()
#4. Functional Annotation----
#Functional Annotations using the likelihood ratio test results 
#This is to measure genes that the overall stage had an effect on as opposed to pairwise comparisons
res_df <- as.data.frame(res_lrt)

# Convert Ensembl IDs to Entrez IDs
# Remove version numbers (e.g., .9 from ENSG00000189221.9)
gene_ids <- rownames(res_df)
gene_ids_clean <- sub("\\..*", "", gene_ids)

# Map to Entrez IDs
gene_map <- bitr(gene_ids_clean, 
                 fromType = "ORF", 
                 toType = c("ENTREZID", "GENENAME"),
                 OrgDb = org.Sc.sgd.db)
# Add the mapping to results
res_df$ORF <- sub("\\..*", "", rownames(res_df))
res_df <- merge(res_df, gene_map, by = "ORF", all.x = TRUE)

# Define significant genes
sig_genes <- res_df %>%
  filter(padj < 0.05 ) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

all_genes <- res_df %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()
#GO ORA analysis
ego_bp <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = F)
head(as.data.frame(ego_bp))
ego_df <- as.data.frame(ego_bp)
dotplot(ego_bp, showCategory = 20, title = "GO Biological Process")
barplot(ego_bp, showCategory = 15, title = "GO Biological Process")
emapplot(pairwise_termsim(ego_bp), showCategory = 30)

#5. Calculating genes upregulated and downregulated based on stage dependence using LRT----
upregulated_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_up <- enrichGO(gene = upregulated_genes,
                      universe = all_genes,
                      OrgDb = org.Sc.sgd.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = F)
downregulated_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < 0 ) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()
ego_bp_down <- enrichGO(gene = upregulated_genes,
                        universe = all_genes,
                        OrgDb = org.Sc.sgd.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = F)

#Showing upregulation and downregulation in relation to the reference(early)
dotplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")

dotplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes")


barplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")

barplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes")

