# PCA_subtype_JCI
R scripts

## Publication
"Tumor subtype defines distinct pathways of molecular and clinical progression in primary prostate cancer"
https://ascopubs.org/doi/full/10.1200/PO.18.00036

## Install
The current R version is 3.6.1. 

Required SVM package:
install.packages("https://cran.r-project.org/src/contrib/Archive/e1071/e1071_1.7-2.tar.gz", repos=NULL)

## Prerequisite

### 1. SPOP mutant signature
Signature is downloaded from Supplementary Table 1a (https://ascopubs.org/doi/suppl/10.1200/PO.18.00036/)

### 2. SPOP signature normalization based on TCGA FPKM expression data
sig212 <- read.table("del2017/SPOP-RNA-classifier/TCGA_333_SPOP_sig_212genes.txt", sep="\t", header=T, check.names=F)

### 3. SPOP mutant status is derived from TCGA PCA study (PMID: 26544944)
Table S1: https://www.cell.com/fulltext/S0092-8674(15)01339-2#supplementaryMaterial
