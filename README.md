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
