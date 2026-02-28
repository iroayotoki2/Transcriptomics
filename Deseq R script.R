library("airway")
library("DESeq2")
library("tximport")
library(ggplot2)

#Load metadata
sampleTable <- read.csv(file.path(dir,"sample_table.csv"),row.names=1)