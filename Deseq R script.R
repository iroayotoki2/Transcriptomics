library("DESeq2")
library("tximport")
library(ggplot2)
library("apeglm")
library("ashr")

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
#Early biofilm as reference
levels(factor(sampleTable$Stage))
#Deseq analysis
dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~Stage)
dds <- DESeq(dds)

resultsNames(dds)
res_mature_vs_early <- results(dds, name="Stage_Mature.Biofilm_vs_Early.biofilm")
res_thin_vs_early   <- results(dds, name="Stage_Thin.Biofilm_vs_Early.biofilm")
res_mature_vs_thin <- results(dds, contrast=c("Stage", "Mature Biofilm", "Thin Biofilm"))

#Likelihood ratio test to see how much stage improve the model
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
#63% of genes significantly affected by stage 
summary(res_lrt)

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
