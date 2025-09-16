# tidy analysis
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

#RNA-seq of TPM
#ovary
CCOC_df <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/CCOC_RNAseq_TPM.csv") # 200

ea_tpm <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/E-MTAB_RNAseq_TPM.csv") # 7

tpm_GSE160692 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE160692_RNA_seq_TPM.csv") # 11

tpm_GSE189553 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE189553_RNAseq_TPM.csv") # 11

tpm_113 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/tpm113_RNAseq_TPM.csv") # 105

tpm_SD1_OC <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/SD1_RNAseq_TPM.csv") # 6

# kidney

tpm_ccRCC <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/ccRCC_RNAseq_TPM.csv")

# Proteomics - already lo2 normalised
#ovary 1

pro_L <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/Proteomic_land_PROT_TPM.csv")

# ES2
prot_es2 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/ES2_PROT_TPM.csv")

# JHOC-5
prot_jhoc5 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/JHOC5_PROT_TPM.csv")

# OVMANA
prot_ovmana <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/OVMANA_PROT_TPM.csv")


#kidney

# OLD -prot_renal <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/RPPA_gene_x_sample_matrix.csv") # too small - old

prot_renal <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/TCGA-KIRC.protein_with_genes2.csv")

#################################
# only get mRNA - checks

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(CCOC_df$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
CCOC_df_mRNA <- CCOC_df %>% filter(Geneid %in% mRNA_genes)

#######
genes <- unique(ea_tpm$Geneid)
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = genes,
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
ea_tpm_mRNA <- ea_tpm %>% filter(Geneid %in% mRNA_genes)

#######

annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_GSE160692$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_GSE160692_mRNA <- tpm_GSE160692 %>% filter(Geneid %in% mRNA_genes)

#######

annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_GSE189553$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_GSE189553_mRNA <- tpm_GSE189553 %>% filter(Geneid %in% mRNA_genes)
#######

annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_113$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_113_mRNA <- tpm_113 %>% filter(Geneid %in% mRNA_genes) 

########
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_SD1_OC$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_SD1_OC_mRNA <- tpm_SD1_OC %>% filter(Geneid %in% mRNA_genes) 

#########
tpm_ccRCC <- tpm_ccRCC %>% rownames_to_column(var = "Geneid")
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_ccRCC$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_ccRCC_mRNA <- tpm_ccRCC %>% filter(Geneid %in% mRNA_genes) 





###############
# RNA-seq only
#############################################
# OVARY
#############################
# sizes
# CCOC_df: 19330   200
# ea_tpm: 38477     7
# tpm_GSE160692: 40405    11
# tpm_GSE189553:56440    11
# tpm_113: 171 105
# tpm_SD1_OC:24427     6

# tpm_ccRCC: 40969   611

tpm_SD1_OC_mRNA <- as.data.frame(tpm_SD1_OC_mRNA)
rownames(tpm_SD1_OC_mRNA) <- make.unique(as.character(tpm_SD1_OC_mRNA$Geneid)) # set row names from colnames
tpm_SD1_OC_mRNA$Geneid <- NULL # delete col

tpm_113_mRNA <- as.data.frame(tpm_113_mRNA)
rownames(tpm_113_mRNA) <- make.unique(as.character(tpm_113_mRNA$Geneid))
tpm_113_mRNA$Geneid <- NULL

tpm_GSE189553_mRNA <- as.data.frame(tpm_GSE189553_mRNA)
rownames(tpm_GSE189553_mRNA) <- make.unique(as.character(tpm_GSE189553_mRNA$Geneid))
tpm_GSE189553_mRNA$Geneid <- NULL

ea_tpm <- as.data.frame(ea_tpm_mRNA)
rownames(ea_tpm_mRNA) <- make.unique(as.character(ea_tpm_mRNA$Geneid))
ea_tpm_mRNA$Geneid <- NULL

tpm_GSE160692_mRNA <- as.data.frame(tpm_GSE160692_mRNA)
rownames(tpm_GSE160692_mRNA) <- make.unique(as.character(tpm_GSE160692_mRNA$Geneid))
tpm_GSE160692_mRNA$Geneid <- NULL

CCOC_df_mRNA <- as.data.frame(CCOC_df_mRNA)
rownames(CCOC_df_mRNA) <- make.unique(as.character(CCOC_df_mRNA$Geneid))
CCOC_df_mRNA$Geneid <- NULL

###############################
# SD1: log transform and convert to data frames and add Geneid column
tpm_SD1_OC_log <- as.data.frame(log2(as.matrix(tpm_SD1_OC_mRNA) + 1))
rownames(tpm_SD1_OC_log) <- rownames(tpm_SD1_OC_mRNA)
tpm_SD1_OC_log$Geneid <- rownames(tpm_SD1_OC_mRNA)

# 113
tpm_113_log <- as.data.frame(log2(as.matrix(tpm_113_mRNA) + 1))
rownames(tpm_113_log) <- rownames(tpm_113_mRNA)
tpm_113_log$Geneid <- rownames(tpm_113_mRNA)

# GSE189553
tpm_GSE189553_log <- as.data.frame(log2(as.matrix(tpm_GSE189553_mRNA) + 1))
rownames(tpm_GSE189553_log) <- rownames(tpm_GSE189553_mRNA)
tpm_GSE189553_log$Geneid <- rownames(tpm_GSE189553_mRNA)

# EA
ea_tpm_log <- as.data.frame(log2(as.matrix(ea_tpm_mRNA) + 1))
rownames(ea_tpm_log) <- rownames(ea_tpm_mRNA)
ea_tpm_log$Geneid <- rownames(ea_tpm_mRNA)

# GSE160692
tpm_GSE160692_log <- as.data.frame(log2(as.matrix(tpm_GSE160692_mRNA) + 1))
rownames(tpm_GSE160692_log) <- rownames(tpm_GSE160692_mRNA)
tpm_GSE160692_log$Geneid <- rownames(tpm_GSE160692_mRNA)

# CCOC
CCOC_df_log <- as.data.frame(log2(as.matrix(CCOC_df_mRNA) + 1))
rownames(CCOC_df_log) <- rownames(CCOC_df_mRNA)
CCOC_df_log$Geneid <- rownames(CCOC_df_mRNA)
#########################

# merge all datasets by Geneid
rna_ovary_merged_log <- merge(tpm_SD1_OC_log, tpm_113_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, tpm_GSE189553_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, ea_tpm_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, tpm_GSE160692_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, CCOC_df_log, by = "Geneid", all = TRUE)

# rownames and clean 
rownames(rna_ovary_merged_log) <- rna_ovary_merged_log$Geneid
rna_ovary_merged_log$Geneid <- NULL
# NAs to0
rna_ovary_merged_log[is.na(rna_ovary_merged_log)] <- 0


######
# variance for each gene
gene_vars <- apply(rna_ovary_merged_log, 1, var)

# top 500 most variable genes
top_genes_500 <- order(gene_vars, decreasing = TRUE)[1:500]

# data to top 500 genes
rna_ovary_top500_log <- rna_ovary_merged_log[top_genes_500, ]

# Z-score 
rna_ovary_top500_z <- t(scale(t(as.matrix(rna_ovary_top500_log))))
rna_ovary_top500_z[is.na(rna_ovary_top500_z)] <- 0

###################
# gene x gene 500
###################
# above 18gb so apply corr after filtering
cor_genes_top <- cor(t(rna_ovary_top500_z), method = "spearman")

Heatmap(rna_ovary_top500_z,
        name = "Spearman",
        show_row_names = FALSE,
        show_column_names = FALSE,
        clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman")),
        clustering_distance_columns = function(x) as.dist(1 - cor(t(x), method = "spearman")),
        column_title = "Top 500 Variable Genes - OCCC RNA-seq",
        row_title = "Genes")


# group by sample name - to see whihc dataset these are from
# group sizes
n_SD1_OC <- ncol(tpm_SD1_OC_log)
n_113 <- ncol(tpm_113_log)
n_GSE189553 <- ncol(tpm_GSE189553_log)
n_ea <- ncol(ea_tpm_log)
n_GSE160692 <- ncol(tpm_GSE160692_log)
n_CCOC <- ncol(CCOC_df_log)

samples <- colnames(rna_ovary_top500_z)

sample_origin <- sapply(samples, function(s) {
  if (s %in% colnames(tpm_SD1_OC_log)) return("SD1_OC")
  else if (s %in% colnames(tpm_113_log)) return("113")
  else if (s %in% colnames(tpm_GSE189553_log)) return("GSE189553")
  else if (s %in% colnames(ea_tpm_log)) return("EA")
  else if (s %in% colnames(tpm_GSE160692_log)) return("GSE160692")
  else if (s %in% colnames(CCOC_df_log)) return("CCOC")
  else return("Unknown")
})
# length(sample_origin) == ncol(rna_ovary_top500_z)


group_colors <- c(
  "SD1_OC" = "#D95F02",
  "113" = "#1B9E77",
  "GSE189553" = "#7570B3",
  "EA" = "#E7298A",
  "GSE160692" = "#66A61E",
  "CCOC" = "#E6AB02"
)

ha <- HeatmapAnnotation(
  Dataset = sample_origin,
  col = list(Dataset = group_colors),
  show_annotation_name = TRUE
)

Heatmap(rna_ovary_top500_z,
        name = "Z-score",
        top_annotation = ha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman")),
        clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman")),
        column_title = "Top 500 Variable Genes - OCCC RNA-seq",
        row_title = "Genes")

range(rna_ovary_top500_z) # -2.086913  6.337701
# strong skew toward positive values
# grey areas: closer to zero: gray/white

##############
#TOP 50
##############

gene_vars <- apply(rna_ovary_merged_log, 1, var)
top_genes_50 <- order(gene_vars, decreasing = TRUE)[1:50]
top_gene_names_50 <- rownames(rna_ovary_merged_log)[top_genes_50]
top_50_z <- rna_ovary_top500_z[top_gene_names_50, ]
# annotation similaryly but for 50 
ha <- HeatmapAnnotation(
  Dataset = sample_origin[1:ncol(top_50_z)],
  col = list(Dataset = group_colors),
  show_annotation_name = TRUE
)

cor_genes_top50 <- cor(t(top_50_z), method = "spearman")
# detect sample and gene clusters
# below to see different study differences - dteect clusters and check the gene expression in diffenret ocnditions, detect outliers in smaples or genes
length(sample_origin) == ncol(top_50_z)  # TRUE
Heatmap(top_50_z, # gene-by-sample expression matrix
        name = "Z-score",
        top_annotation = ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman")),
        clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman")),
        column_title = "Top 50 Variable Genes - OCCC RNA-seq",
        row_title = "Genes")

# for gene to gene relationships
Heatmap(cor_genes_top50, # gene-by-gene correlation matrix
        name = "Spearman", 
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = function(x) as.dist(1 - x),
        clustering_distance_columns = function(x) as.dist(1 - x),
        column_title = "Top 50 Genes Correlation (OCCC)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_genes_top50[i, j]), x, y, gp = gpar(fontsize = 6))
        })


###################
# RNA-seq only
##########################
# KIDNEY
#########################

tpm_ccRCC <- as.data.frame(tpm_ccRCC)
tpm_ccRCC <- tpm_ccRCC[!is.na(tpm_ccRCC$Geneid), ]  # remove NA Geneids
rownames(tpm_ccRCC) <- make.unique(as.character(tpm_ccRCC$Geneid))
tpm_ccRCC$Geneid <- NULL  # remove column after using as rownames

# logtr
tpm_ccRCC_log <- log2(data.matrix(tpm_ccRCC) + 1)

# rank genes by variance **before** z-scoring
gene_vars_kidney <- apply(tpm_ccRCC_log, 1, var)
top_genes_kidney <- order(gene_vars_kidney, decreasing = TRUE)[1:50]

#  longtr data to top variable genes
rna_kidney_log_top <- tpm_ccRCC_log[top_genes_kidney, ]

# zscore normalise (gene)
rna_kidney_z_top <- t(scale(t(rna_kidney_log_top)))
rna_kidney_z_top[is.na(rna_kidney_z_top)] <- 0  # optional, remove NA from low variance

cor_genes_kidney <- cor(t(rna_kidney_z_top), method = "spearman")

# rna_kidney_z_top if this is applied to heatmap: Which genes are up/downregulated in each sample
# cor_genes_kidney: Co-expression between genes: genes behave similarly → potential for same biological pathways

Heatmap(cor_genes_kidney,
        name = "Spearman Correlation",
        show_row_names = TRUE,
        show_column_names = TRUE,
        clustering_distance_rows = function(x) as.dist(1 - x),
        clustering_distance_columns = function(x) as.dist(1 - x),
        column_title = "Top 50 Variable Genes Correlation (ccRCC)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_genes_kidney[i, j]), x, y, gp = gpar(fontsize = 6))
        })


# I can also look into bigger span rather than doing heatmap just get common genes 
# and then find biomarkers like top 500

##############
###########################
# DE fro log2tr values (before zscore)
###########################
##############

rna_ovary_merged_log

rownames(rna_ovary_merged_log) <- rna_ovary_merged_log$Geneid
rna_ovary_merged_log$Geneid <- NULL
rna_ovary_merged_log[is.na(rna_ovary_merged_log)] <- 0
rna_ovary_merged_log[] <- lapply(rna_ovary_merged_log, as.numeric)
# remove genes with zero variance
rna_ovary_filtered <- rna_ovary_merged_log[apply(rna_ovary_merged_log, 1, var) != 0, ]

pca_result <- prcomp(t(rna_ovary_filtered), scale. = TRUE)

pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")

head(pca_result$x)
head(pca_result$rotation)
summary(pca_result)
# PC1 - 27.16%
# PC2 - 16.60%
# PC3 - 12.44%
# healthy distribution of variance, not noise-dominated data or weird structure
# PC1–PC3 explain ~56% of the variance - good - biologically meaningful variation not just noise
# PCA picked up major patterns, which could be disease subtype, treatment, sample source, or batch so my data is informative
plot(pca_result, type = "l", main = "Scree Plot")

# only Ovary
# expression matrix
exprs <- rna_ovary_merged_log
# sample metadata
samples <- colnames(exprs)

# assign group name
#colnames(tpm_SD1_OC_log) <- paste0("SD1_", colnames(tpm_SD1_OC_log))
colnames(tpm_GSE160692_log) <- paste0("GSE160692_", colnames(tpm_GSE160692_log))


ea_samples <- c("ES-2", "JHOC-5", "OVISE", "OVMANA", "OVTOKO", "RMG-I", "TOV-21G")

group <- ifelse(samples %in% ea_samples, "EA",
                ifelse(samples %in% c("C1", "C2", "C3", "C4", "C5", "C6"), "SD1",
                       ifelse(grepl("- RPKM$", samples), "S113",
                              ifelse(grepl("^CCC_", samples), "GSE189553",
                                     ifelse(grepl("^OVA[0-9]+", samples), "GSE160692",
                                            ifelse(grepl("IGO_07456", samples), "CCOC", "Unknown"))))))


#  metadata data.frame
sample_info <- data.frame(Sample = samples, Group = factor(group))
# design matrix
design <- model.matrix(~0 + sample_info$Group)
colnames(design) <- levels(sample_info$Group)

# filtering low expressed genes
# tpm already normalised so, can skip if not or do:
keep <- rowMeans(exprs) > 1
exprs_filtered <- exprs[keep, ]

# apply limma
fit <- lmFit(exprs_filtered, design)

# add in all
contrast_matrix <- makeContrasts(
  #SD1
  CCOC_vs_SD1 = CCOC - SD1,
  EA_vs_SD1 = EA - SD1,
  GSE160692_vs_SD1 = GSE160692 - SD1,
  GSE189553_vs_SD1 = GSE189553 - SD1,
  S113_vs_SD1 = S113 - SD1,
  #113
  S113_vs_GSE189553 = S113 - GSE189553,
  S113_vs_GSE160692 = S113 - GSE160692,
  S113_vs_EA = S113 - EA,
  S113_vs_CCOC = S113 - CCOC,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results_ccoc <- topTable(fit2, coef = "CCOC_vs_SD1", number = Inf)
results_ea <- topTable(fit2, coef = "EA_vs_SD1", number = Inf)
results_gse160692 <- topTable(fit2, coef = "GSE160692_vs_SD1", number = Inf)


head(results_ccoc) 
head(results_ea)
head(results_gse160692)

##########################
#merged kidney and ovary DE
##########################
#rna_ovary_merged_log vs tpm_ccRCC_log

# convert both to data frames
rna_ovary_merged_log <- as.data.frame(rna_ovary_merged_log)
tpm_ccRCC_log <- as.data.frame(tpm_ccRCC_log)

# add Geneid as a column from rownames
rna_ovary_merged_log$Geneid <- rownames(rna_ovary_merged_log)
tpm_ccRCC_log$Geneid <- rownames(tpm_ccRCC_log)

# merge on Geneid
merged_data <- merge(rna_ovary_merged_log, tpm_ccRCC_log, by = "Geneid", all = FALSE)

# move Geneid back to rownames and remove the column
rownames(merged_data) <- merged_data$Geneid
merged_data$Geneid <- NULL

##### wihtout limma norm 
# convert to expression matrix
#exprs <- as.matrix(merged_data)
# checks of str
#str(exprs)
#head(exprs[, 1:5])

#samples <- colnames(exprs)
# name the groups 
#group <- ifelse(grepl("^TCGA-", samples), "ccRCC", "OCCC")

#sample_info <- data.frame(Sample = samples, Group = factor(group))
#design <- model.matrix(~0 + sample_info$Group)
#colnames(design) <- levels(sample_info$Group)
# filter low
#keep <- rowMeans(exprs) > 1
#exprs_filtered <- exprs[keep, ]
############################

################################
# try apply limma norm
##sample info and design 
samples <- colnames(exprs)

# Label groups from sample IDs (adjust the grep rule if needed)
group <- ifelse(grepl("^TCGA-", samples), "ccRCC", "OCCC")

# Lock the level order so the contrast direction is clear
sample_info <- data.frame(Sample = samples,
                          Group  = factor(group, levels = c("OCCC", "ccRCC")))

# Ensure columns of exprs match sample_info order (defensive)
exprs <- exprs[, sample_info$Sample, drop = FALSE]

design <- model.matrix(~ 0 + Group, data = sample_info)
colnames(design) <- levels(sample_info$Group)  # "OCCC", "ccRCC"

## filter low expression (TPM ≥ 1 in ≥ 20% of samples) 
# Back-transform log2(TPM+offset) -> TPM.
# If your log used +1 as offset, this is correct. If you used a different offset, change it below.
TPM <- 2^exprs - 1
keep <- rowMeans(TPM >= 1) >= 0.20
exprs <- exprs[keep, , drop = FALSE]

##Optional between-sample normalization (quantile)
exprs_norm <- normalizeBetweenArrays(exprs, method = "quantile")

## fit limma model on log-TPM intensities 
fit  <- lmFit(exprs_norm, design)

# Contrast: ccRCC vs OCCC (positive logFC means higher in ccRCC)
contr <- makeContrasts(ccRCC_vs_OCCC = ccRCC - OCCC, levels = design)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

res  <- topTable(fit2, coef = "ccRCC_vs_OCCC", number = Inf, adjust.method = "BH")

# (Optional) quick QC plots
plotDensities(exprs, main = "Before quantile")
plotDensities(exprs_norm, main = "After quantile")
plotSA(fit2)
plotMD(fit2, coef = 1)

#################################

# VIEWING LIMMA NORM RESULTS 


fc_thresh <- 1
fdr_thresh <- 0.05

res$Significance <- "Not Sig"
res$Significance[res$adj.P.Val < fdr_thresh & res$logFC > fc_thresh] <- "Up"
res$Significance[res$adj.P.Val < fdr_thresh & res$logFC < -fc_thresh] <- "Down"

ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot: ccRCC vs OCCC") +
  theme_minimal()
#inspect geens
head(res, 20)                           # top 20
sig <- subset(res, adj.P.Val < 0.05)    # significant set
write.csv(sig, "DE_genes_qnorm.csv", row.names = FALSE)

library(ComplexHeatmap)
library(circlize)

## --- pick significant genes from limma results ------------------------------
# adjust threshold / topN as you like
sig_genes <- subset(res, adj.P.Val < 0.05)
sig_genes <- sig_genes[order(sig_genes$adj.P.Val, -abs(sig_genes$logFC)), , drop = FALSE]

# cap the heatmap to a manageable number of rows if needed
topN <- min(100, nrow(sig_genes))  # show up to 100 genes
sig_gene_names <- rownames(sig_genes)[seq_len(topN)]

## build exp matrix
# use your normalized log-TPM matrix
mat <- exprs_norm

# keep only genes present in the matrix
sig_gene_names <- intersect(sig_gene_names, rownames(mat))
heatmap_data <- mat[sig_gene_names, , drop = FALSE]

# ensure columns (samples) align with sample_info
heatmap_data <- heatmap_data[, sample_info$Sample, drop = FALSE]

## (optional) row scaling to highlight relative changes 
#z-score each gene across samples
z <- t(scale(t(heatmap_data)))
# remove any rows that became NA due to zero variance
z <- z[complete.cases(z), , drop = FALSE]

## anno
grp <- droplevels(sample_info$Group)

# color map for groups
grp_levels <- levels(grp)
grp_cols <- setNames(hcl.colors(length(grp_levels), "Dark2"), grp_levels)

col_ha <- HeatmapAnnotation(
  Group = grp,
  col = list(Group = grp_cols),
  annotation_name_side = "left"
)


col_fun <- colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))

Heatmap(
  z,
  name = "z-score (log2 TPM)",
  top_annotation = col_ha,
  show_row_names = nrow(z) <= 50,        # avoid clutter if many genes
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "Significant DE genes",
  heatmap_legend_param = list(title = "z-score"),
  use_raster = TRUE, raster_quality = 2
)

# gene x gene
library(ComplexHeatmap)
library(circlize)

# order by adjusted p-value and effect size
sig_genes <- res[order(res$adj.P.Val, -abs(res$logFC)), ]
top50 <- head(rownames(sig_genes), 50)

mat <- exprs_norm[top50, , drop = FALSE]

cor_mat <- cor(t(mat), method = "spearman")  # gene × gene correlations

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(
  cor_mat,
  name = "Spearman", # spearman r
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_title = "Top 50 DE Genes: ccRCC vs OCCC", # gene x gene
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 6))
  }
)


# RANKED LIST
#################
# GSEA - Ranked list of all genes - logFC
# ranked gene list for GSEA (NO filtering)
gene_list <- res$logFC
names(gene_list) <- rownames(res)

# Drop NAs and sort
gene_list <- sort(na.omit(gene_list), decreasing = TRUE)

# Map SYMBOL to ENTREZ IDs
gene_df <- bitr(names(gene_list), fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db) # 0.11% of input gene IDs are fail to map

# Align and rename
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID

# Run ranked GSEA
gsea_result <- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

head(gsea_result)

write.table(data.frame(SYMBOL = names(gene_list), logFC = gene_list),
            file = "ranked_gene_listGSEA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)


#---------
# keep pathwya enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

head(kegg_enrich)

# Bar plot of top GO terms
barplot(go_enrich, showCategory = 20)

# Dot plot
dotplot(go_enrich, showCategory = 20)

# KEGG plot
barplot(kegg_enrich, showCategory = 20)



###########

# apply limma
fit <- lmFit(exprs_filtered, design)

# ccRCC vs OCCC
contrast_matrix <- makeContrasts(ccRCC - OCCC, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf, adjust.method = "fdr")
head(results)
#all
write.csv(results, file = "ccRCC_vs_OCCC_DE_results.csv")


# significant genes # "results" for volcano too
sig_genes <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(sig_genes) # 1804 numebr of sig genes

# --OLD Volcano--assess DE gene significance and effect size 
with(results, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano Plot",
                   xlab="log2 Fold Change", ylab="-log10 FDR"))
abline(h = -log10(0.05), col = "red", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)
####

# top 50 significant genes
top_genes <- rownames(sig_genes)[1:50]
heatmap_data <- exprs_filtered[top_genes, ]

# Scale and plot
pheatmap(heatmap_data, scale = "row", show_rownames = TRUE, show_colnames = FALSE,
         main = "Top 50 DE Genes (ccRCC vs OCCC)")


# another heatmap with gene names
sig_gene_names <- rownames(sig_genes)
# expression matrix
heatmap_data <- exprs_filtered[sig_gene_names, ]
# sample info (same order as exprs)
group_factor <- factor(group)
# annotation for columns
col_ha <- HeatmapAnnotation(
  Group = group_factor,
  col = list(Group = c("ccRCC" = "red", "Ovary" = "blue"))
)

Heatmap(
  heatmap_data,
  name = "log2 TPM",
  top_annotation = col_ha,
  show_row_names = nrow(heatmap_data) <= 50,  # avoid clutter if too many genes
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "Significant DE Genes: ccRCC vs Ovary",
  heatmap_legend_param = list(title = "Expression")
)

# select top 50
top_n <- 50
sig_genes_top <- head(sig_genes[order(sig_genes$adj.P.Val), ], top_n)
top_gene_names <- rownames(sig_genes_top)


heatmap_data <- exprs_filtered[top_gene_names, ]
# already defined group factro and col_ha so no need to write here again

pdf("Top50_DEG_Heatmap.pdf", width = 8, height = 10)
Heatmap(
  heatmap_data,
  name = "log2 TPM",
  top_annotation = col_ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "Top 50 DE Genes: ccRCC vs Ovary",
  heatmap_legend_param = list(title = "Expression")
)
dev.off()

# prev haev too many sample
#VOLCANO
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'ccRCC vs Ovary',
                subtitle = "Volcano Plot",
                legendLabels = c("NS", "Log2FC", "Adj.P", "Both"),
                legendPosition = 'right'
)


# NEW VOLCANO W SIG LABEL
# thresholds
fc_thresh <- 1
fdr_thresh <- 0.05

results$Significance <- "Not Sig"
results$Significance[results$adj.P.Val < fdr_thresh & results$logFC > fc_thresh] <- "Up"
results$Significance[results$adj.P.Val < fdr_thresh & results$logFC < -fc_thresh] <- "Down"

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot: ccRCC vs OCCC") +
  theme_minimal()





################
##########################

# GSEA
######################
library(clusterProfiler)
library(org.Hs.eg.db)

# Over-Representation Analysis- Filtered list of significant genes
sig_genes <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
gene_symbols <- rownames(sig_genes)

# map gene symbols to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db) # some didnt mapped - 0.1%

entrez_ids <- gene_df$ENTREZID

go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",        # BP = Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

head(go_enrich)


write.csv(as.data.frame(go_enrich), "GO_enrichment_ccRCC_vs_OCCC.csv", row.names = FALSE)
#################
# RANKED LIST
#################
# GSEA - Ranked list of all genes - logFC
# ranked gene list for GSEA (NO filtering)
gene_list <- results$logFC
names(gene_list) <- rownames(results)

# Drop NAs and sort
gene_list <- sort(na.omit(gene_list), decreasing = TRUE)

# Map SYMBOL to ENTREZ IDs
gene_df <- bitr(names(gene_list), fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db) # 0.11% of input gene IDs are fail to map

# Align and rename
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID

# Run ranked GSEA
gsea_result <- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

head(gsea_result)

write.table(data.frame(SYMBOL = names(gene_list), logFC = gene_list),
            file = "ranked_gene_listGSEA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)


#---------
# keep pathwya enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

head(kegg_enrich)

# Bar plot of top GO terms
barplot(go_enrich, showCategory = 20)

# Dot plot
dotplot(go_enrich, showCategory = 20)

# KEGG plot
barplot(kegg_enrich, showCategory = 20)

###########
# I cna also use ranked gene list # ranks by logFC #
gene_list <- results$logFC
names(gene_list) <- rownames(results)

# convert to Entrez IDs
gene_df_full <- bitr(names(gene_list), fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

gene_list <- gene_list[gene_df_full$SYMBOL]
names(gene_list) <- gene_df_full$ENTREZID

# Sort
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_result <- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     verbose = FALSE)

dotplot(gsea_result, showCategory = 20)
###
#2nd ranked list # ranks by t-statistic #
# assuming 'results' has rownames = SYMBOL and columns: t, logFC, adj.P.Val
results$SYMBOL <- rownames(results)
rank_map <- bitr(results$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ranks_df <- merge(rank_map, results[,c("SYMBOL","t")], by="SYMBOL", all.x=TRUE)
ranks <- ranks_df$t; names(ranks) <- ranks_df$ENTREZID
ranks <- sort(ranks, decreasing=TRUE)

gsea_bp <- gseGO(geneList=ranks, OrgDb=org.Hs.eg.db, ont="BP", verbose=FALSE)
dotplot(gsea_bp, showCategory = 20)
# for a signal-to-noise–type score more robust than the other one


#################
# finding oncogenes
#################

#  oncogene list 
oncokb <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/cancerGeneList.tsv", show_col_types = FALSE)

# clean col
oncokb_clean <- oncokb[, c("Hugo Symbol", "Gene Type", "OncoKB Annotated")]
colnames(oncokb_clean)[1] <- "Gene"  # rename 

# adding gene names to the sig_genes dataframe as a column
sig_genes$Gene <- rownames(sig_genes)
# merge DE results 
sig_annotated <- merge(sig_genes, oncokb_clean, by = "Gene", all.x = TRUE)

# Check result
head(sig_annotated)
# Only keep DE genes that are oncogenes or tumor suppressors
sig_onco_tsg <- sig_annotated[!is.na(sig_annotated$`Gene Type`), ]

#sig_oncogenes <- subset(sig_onco_tsg, `Gene Type` == "ONCOGENE")
#sig_tsgs <- subset(sig_onco_tsg, `Gene Type` == "TSG")
write.csv(sig_onco_tsg, "significant_oncogenes_tsgs.csv", row.names = FALSE)


ggplot(sig_onco_tsg, aes(x = logFC, y = -log10(adj.P.Val), color = `Gene Type`)) +
  geom_point() +
  geom_text(aes(label = Gene), vjust = 1, hjust = 1, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "DE Oncogenes & TSGs", x = "log2 Fold Change", y = "-log10(FDR)")



# significant oncogenes clusterprofiler
oncogenes <- oncokb %>%
  filter(`Gene Type` == "ONCOGENE") %>%
  pull(`Hugo Symbol`) %>%
  toupper()
# to get sig_genes_filtered run kinases first
sig_oncogenes <- sig_genes_filtered %>% filter(Gene %in% oncogenes)
# Oncogenes
oncogene_entrez <- bitr(sig_oncogenes$Gene, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Oncogenes
ego_onco <- enrichGO(gene = oncogene_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = TRUE)

# Oncogenes kegg
ekegg_onco <- enrichKEGG(gene = oncogene_entrez$ENTREZID,
                         organism = 'hsa', pvalueCutoff = 0.05)

barplot(ego_onco, showCategory = 10, title = "GO Enrichment - Oncogenes")
dotplot(ego_onco, showCategory = 10, title = "GO Dotplot - Oncogenes")
dotplot(ekegg_onco, showCategory = 10, title = "KEGG - Oncogenes")

cnetplot(ego_onco, categorySize = "pvalue")

###############
# finding kinases
###############


# Example list of human kinases (a small subset for demo)
kinases <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/kinase_basic.xlsx", col_names = TRUE)
colnames(kinases)

top_genes <- head(results[order(results$adj.P.Val), ], 1000) # use larger
# mostly pseudogenens adn non codinf RNAs so cannot find kinases
# Your significant genes
sig_kinases <- kinases[kinases$Offical_gene_symbol %in% top_genes, ]
head(sig_kinases)


# reorder and filter significant genes with logFC threshold
sig_genes_filtered <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
sig_genes_filtered$Gene <- toupper(rownames(sig_genes_filtered))

# match against kinase genes
kinases$Offical_gene_symbol <- toupper(kinases$Offical_gene_symbol)

sig_kinases <- sig_genes_filtered[sig_genes_filtered$Gene %in% kinases$Offical_gene_symbol, ]
nrow(sig_kinases)
head(sig_kinases)

write.csv(sig_kinases, "significant_kinases.csv", row.names = FALSE)

# expression matrix for  kinases
kinase_expr <- exprs_filtered[rownames(exprs_filtered) %in% sig_kinases$Gene, ]
# annotation
ha <- HeatmapAnnotation(Group = sample_info$Group)
# Plot
Heatmap(kinase_expr, top_annotation = ha, show_row_names = TRUE, show_column_names = FALSE,
        cluster_rows = TRUE, cluster_columns = TRUE,heatmap_legend_param = list(title = "Expression"))


#####clusterprofiler#####
# significant kinases
kinase_symbols <- kinases$Offical_gene_symbol

sig_kinases <- sig_genes_filtered %>%
  filter(Gene %in% kinase_symbols)

# Kinases
kinase_entrez <- bitr(sig_kinases$Gene, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Kinases
ego_kinase <- enrichGO(gene = kinase_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, readable = TRUE)


# Kinases kegg
ekegg_kinase <- enrichKEGG(gene = kinase_entrez$ENTREZID,
                           organism = 'hsa', pvalueCutoff = 0.05)

barplot(ego_kinase, showCategory = 10, title = "GO Enrichment - Kinases")
dotplot(ego_kinase, showCategory = 10, title = "GO Dotplot - Kinases")
dotplot(ekegg_kinase, showCategory = 10, title = "KEGG - Kinases")

cnetplot(ego_kinase, categorySize = "pvalue")


write.csv(as.data.frame(ego_onco), "GO_oncogene_enrichment.csv")
write.csv(as.data.frame(ego_kinase), "GO_kinase_enrichment.csv")

##############################################
# REACTOME - RNA O vs K
######################
# install packages once if needed:
#BiocManager::install("ReactomePA")

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Suppose your limma/DE table is in `res` with columns: SYMBOL, logFC, adj.P.Val, t, P.Value
# If symbols are rownames, do: res$SYMBOL <- rownames(res)

# significant genes for ORA
sig <- results %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% distinct(SYMBOL)

# map SYMBOL → Entrez (required by ReactomePA)
entrez_sig <- bitr(sig$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

# 3) build a ranked vector for GSEA (names = Entrez, values = ranking metric)
ranks_df <- bitr(res$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>%
  inner_join(res, by = "SYMBOL") %>%
  select(ENTREZID, t) # you can use t or logFC

ranks <- ranks_df$t
names(ranks) <- ranks_df$ENTREZID
ranks <- sort(ranks, decreasing = TRUE)

#########
# Reactome ORA (over-representation) + plots

# sig genes
sig_syms <- results$SYMBOL[results$adj.P.Val < 0.05 & abs(results$logFC) > 1]

# map to Entrez
sig_entrez <- bitr(sig_syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

# set universe = all tested genes
bg_entrez <- bitr(results$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

ora <- enrichPathway(
  gene         = unique(sig_entrez),
  organism     = "human",
  universe     = unique(bg_entrez),
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable     = TRUE
)


dotplot(ora, showCategory = 20, title = "Reactome ORA")
barplot(ora, showCategory = 20, title = "Reactome ORA")

# network-style plots
ora_sim <- pairwise_termsim(ora)
emapplot(ora_sim, showCategory = 30, layout = "kk") # 
cnetplot(ora, showCategory = 10, categorySize = "pvalue", # 11 unlabeled data points 
         foldChange = setNames(results$logFC, results$SYMBOL))
# 17 unlabeled data points (too many overlaps)

###############

# ranked vector using limma t-stat (recommended over logFC)
rank_map <- bitr(results$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ranks_df <- inner_join(rank_map, results[, c("SYMBOL","t")], by="SYMBOL") %>% 
  filter(!is.na(t)) %>% 
  group_by(ENTREZID) %>%  # deduplicate Entrez, keep strongest |t|
  slice_max(order_by = abs(t), n = 1) %>% 
  ungroup()

ranks <- ranks_df$t
names(ranks) <- ranks_df$ENTREZID
ranks <- sort(ranks, decreasing = TRUE)

# run GSEA (Reactome)
gsea <- gsePathway(
  geneList     = ranks,
  organism     = "human",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)


dotplot(gsea, showCategory = 20, split = ".sign", title = "Reactome GSEA") # BP = Gene Ontology Biological Process
# split by NES sign - active vs suppressed
dotplot(gsea, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) # group the results by the variable .sign and split panel

ridgeplot(gsea, showCategory = 20, fill = "pvalue")

library(enrichplot) 
# Enrichment plot for the most significant term
gdf <- as.data.frame(gsea)
top_id <- gdf$ID[which.min(gdf$p.adjust)]
gseaplot2(gsea, geneSetID = top_id, title = gdf$Description[gdf$ID == top_id])
# pathway is positively enriched (NES > 0)
# covers ribosome/ER targeting via SRP—i.e., increased protein synthesis/secretory pathway engagement
# plot shows that SRP-dependent cotranslational ER targeting is significantly upregulated in the group that defines positive ranks likely ccRCC 








##################################################
##################################################
# PROTEOMICS only
################

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


pro_L
pro_L <- as.data.frame(pro_L)
rownames(pro_L) <- make.unique(as.character(pro_L$Geneid))
pro_L$Geneid <- NULL

pro_L_log <- as.data.frame(log2(data.matrix(pro_L) + 1))

pro_L_log$Geneid <- rownames(pro_L_log)

# Z-score 
rna_ovary_top500_z <- t(scale(t(as.matrix(rna_ovary_top500_log))))
rna_ovary_top500_z[is.na(rna_ovary_top500_z)] <- 0


# merge all datasets by Geneid of prot
rna_ovary_merged_log <- merge(tpm_SD1_OC_log, tpm_113_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, tpm_GSE189553_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, ea_tpm_log, by = "Geneid", all = TRUE)
rna_ovary_merged_log <- merge(rna_ovary_merged_log, tpm_GSE160692_log, by = "Geneid", all = TRUE)










