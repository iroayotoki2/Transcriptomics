# Differential Gene Expression of Yeast Biofilm during Wine Aging

Differential gene expression is the process by which cells with the same genome activate different genes for different adaptations and functions, and studying it allows us to better understand cell differentiation, gene regulation, and adaptation to environments (Gilbert, 2000). The purpose of this project is to research the best methods, software, and parameters for gene expression analysis and to use bulk transcriptomics data to study differential gene expression in yeast biofilm during wine aging.

Velum refers to a surface film formed by yeast (*Saccharomyces cerevisiae*) that grow on top of a liquid, especially in winemaking. It is important to study because it dictates the release of aromatic compounds, shapes the wine’s mouthfeel, and defines its stability and longevity. It also shows how the yeast adapts to high ethanol environments, consumes ethanol and glycerol, and forms a protective, metabolically active biofilm. This helps to optimize the process of producing quality wine and to understand the biological mechanisms involved (Mardanov et al., 2020).

In order to properly study how gene expression affects this process, we must use appropriate and effective tools and software for the data provided to achieve accurate results.

## Methods Comparison

A 2021 study carried out by NASA GeneLab shows FastQC as a reliable quality control tool for short-read sequencing data, as it provides information that allows users to assess sample and sequencing quality, including base statistics, per-base sequencing quality, per-sequence quality scores, per-base sequence content, per-base GC content, per-sequence GC content, per-base N content, sequence length distributions, sequence duplication levels, overrepresented sequences, and k-mer content (Overbey et al., 2021). The outputs from this tool can also be combined to produce a multiqc report comapring structure of multple sequence files. 

Another 2024 study showed STAR to be a superior aligner with an accuracy of 90% when compared to other aligners for RNA-Seq data (Coxe et al., 2024).

A comparison study of quantification tools showed RSEM to be the best quantifier when compared with tools like Salmon and NanoString, with low absolute mean error values (Germain et al., 2016).

Finally, an evaluation of RNA-seq differential analysis methods showed that DESeq2 performs best when the sample size is six or higher per group, offering the best balance of FDR control, power, and stability (Li et al., 2022).

This is optimal for our data structure, as we have nine different samples, and these methods can provide strong statistical results.
# Methods

## Data Acquisition

The dataset used for this analysis was obtained from the NCBI database under BioProject accession PRJNA592304. It consisted of nine samples in SRR file format representing the three stages of biofilm development (Early, Thin, and Mature biofilm). SRR files were converted to FASTQ files using SRA Toolkit Release 3.3.0 (The Sequence Read Archive (SRA), n.d.). The reference genome and annotation for *Saccharomyces cerevisiae* were obtained from the NCBI FTP website in FASTA and GTF format.

## Quality Control

Quality control checks were performed using FastQC v0.12.1 for each individual FASTQ file and were viewed collectively using MultiQC v1.33 (Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data, n.d.; Ewels et al., 2016).

## Alignment

Alignment was performed using STAR v2.7.11b and involved two steps. The first step was generating the genome index using `--runMode genomeGenerate`. The second step involved aligning reads, with the key option being `--quantMode TranscriptomeSAM` to produce BAM files aligned to the transcriptome (Dobin et al., 2013).

## Quantification

Quantification was performed using RSEM v1.3.3 by preparing the reference with the `rsem-prepare-reference` command using the reference FASTA and GTF files. This reference was then supplied to the `rsem-calculate-expression` command to quantify gene counts for the aligned files, using the `--bam` option to specify input. This produced gene-level and isoform-level count files (Li & Dewey, 2011).

All bash shell scripting code can be found in the [transcriptomics.sh](transcriptomics.sh) script file.

## Differential Analysis

Differential expression analysis was performed using DESeq2 v1.44.0 in R v4.5.1. Gene count files generated from quantification were imported using tximport v3.21. DESeq2 was used to produce pairwise comparisons among the three stages (Early, Thin, and Mature). The model was also subjected to a Likelihood Ratio Test (LRT) to assess how inclusion of Stage improved model fit (Love et al., 2014).

## Functional Annotation

Functional annotation was performed using the clusterProfiler R package v4.14.3 to conduct Gene Ontology (GO) overrepresentation analysis. The results data frame from the differential analysis LRT was used to identify biological processes significantly associated with Stage and to compute upregulated and downregulated genes (Wu et al., 2021).

## Data Visualization

Data structure and results were visualized using R packages including ggplot2 v4.0.2, pheatmap v1.0.13, and enrichplot v1.31.4 (Bioconductor - Enrichplot, n.d.; “The Grammar of Graphics,” 2005; Kolde, 2025).
All R Code can be found in the [Deseq R Script](Deseq%20R%20script.R)
