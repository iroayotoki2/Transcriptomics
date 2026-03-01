library("DESeq2")
library("tximport")
library(ggplot2)

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
