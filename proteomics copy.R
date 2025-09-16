##################################################
# PROTEOMICS only
################
library(tidyverse)
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(limma)
library(readxl)
library(EnhancedVolcano)
library(biomaRt)

# Proteomics - already lo2 normalised
#ovary 1

pro_L <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/Proteomic_land_PROT.csv") # no norm

# ES2
prot_es2 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/ES2_PROT.csv") #norm

# JHOC-5
prot_jhoc5 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/JHOC5_PROT.csv") #norm 

# OVMANA
prot_ovmana <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/OVMANA_PROT.csv") #norm


#kidney

# OLD -prot_renal <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/RPPA_gene_x_sample_matrix.csv") # too small - old

prot_renal <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/TCGA-KIRC.protein_with_genes2.csv") 



# OVARY
# pro_L : 40 genes, 5 samples - need log2 and z score
#prot_es2: 144 genes, 900 samples -already in z score - so cannot apply log2tr
#prot_jhoc5 : 144 genes 900 samples - already zscore norm
#prot_ovmana: 144 genes 900 smaples - already zscore nrom

# Kidney
# prot_renal: 430 genes 57 samples - already zscore norm


# comparing correlation heatmaps or variance-based rankings — 
# smaller sample sizes lead to less stable correlations and variance estimates

# doesn’t make sense to use variance for filtering after z-score, 
# because z-scoring normalizes each gene to have variance = 1
#1. Log2-transform
#2. Select top variable genes: after log transformation, it is generally safe to calculate and rank variance per gene — even if sample sizes are imbalanced
#3. Then z-score



# Z-score
pro_L_z <- t(scale(t(as.matrix(pro_L))))
pro_L_z <- as.data.frame(pro_L_z)

pro_L_z[is.na(pro_L_z)] <- 0

pro_L_z$Geneid <- rownames(pro_L_z)

# (optional: reorder so Geneid is the first column)
pro_L_z <- pro_L_z[, c("Geneid", setdiff(names(pro_L_z), "Geneid"))]




pro_L <- as.data.frame(pro_L)
rownames(pro_L) <- make.unique(as.character(pro_L$Geneid))
pro_L$Geneid <- NULL

# Z-score 
pro_L_z <- t(scale(t(as.matrix(pro_L))))
pro_L_z[is.na(pro_L)] <- 0
pro_L_z <- as.data.frame(pro_L_z)

# merge all datasets by Geneid of prot
pro_L_merged <- merge(prot_es2, pro_L_z, by = "Geneid", all = TRUE)
pro_L_merged <- merge(pro_L_merged, prot_es2, by = "Geneid", all = TRUE)
pro_L_merged <- merge(pro_L_merged, prot_jhoc5, by = "Geneid", all = TRUE)
pro_L_merged <- merge(pro_L_merged, prot_ovmana, by = "Geneid", all = TRUE)

# why ARHGEF40 is NA














