# updated dataset and analysis script 

library(tidyverse)
library(readxl)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("GenomicRanges", force = TRUE)
BiocManager::install("SummarizedExperiment", force = TRUE)
#BiocManager::install("biomaRt")
#library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(grid)
library(ggplot2)
library(ggrepel)
library(tibble)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)


#####################################
# Molecular Subclasses of Clear Cell Ovarian Carcinoma and Their Impact on Disease Behavior and Outcomes
# (1) sample 211 primary - tpm
#####################################

CCOC_rds_file <- readRDS("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/CCOC.RNAseq.normalized.rds")

# remove _PAR_Y
head(rownames(CCOC_rds_file))
rownames(CCOC_rds_file) <- sub("_PAR_Y$", "", rownames(CCOC_rds_file))
# remove verison
rds_version_ids <- rownames(CCOC_rds_file)
rds_gene_ids <- sub("\\..*$", "", rds_version_ids)
# all are ensg -> gene symbol (onyl mrna)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

bm_results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id", values = rds_gene_ids, mart = ensembl)

bm_results_mrna <- bm_results %>% filter(gene_biotype == "protein_coding")
bm_results_mrna <- bm_results_mrna %>% filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

# merge back
mapping_df <- data.frame(ensembl_with_version = rds_version_ids, ensembl_gene_ids = rds_gene_ids)
# only protein-coding genes to mapping_df
colnames(mapping_df)[colnames(mapping_df) == "ensembl_gene_ids"] <- "ensembl_gene_id"
mapping_df <- merge(mapping_df, bm_results_mrna, by = "ensembl_gene_id",all.x = FALSE)

CCOC_rds_file <- CCOC_rds_file[mapping_df$ensembl_with_version, ]

CCOC_df <- as.data.frame(CCOC_rds_file)
# duplicates
#dup_genes <- duplicated(mapping_df$hgnc_symbol)

CCOC_df$hgnc_symbol <- mapping_df$hgnc_symbol
CCOC_df <- CCOC_df[, c(ncol(CCOC_df), 1:(ncol(CCOC_df)-1))]
# take mean of duplicates- mean for each numeric column separately
CCOC_df <- CCOC_df %>% group_by(hgnc_symbol) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% ungroup()

CCOC_df <- as.data.frame(CCOC_df)
# nrow = 19,330
# chaneg col name
colnames(CCOC_df) <- trimws(colnames(CCOC_df))
colnames(CCOC_df)[colnames(CCOC_df) == "hgnc_symbol"] <- "Geneid"

write.csv(CCOC_df,
          file = "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/CCOC_RNAseq_TPM.csv",
          row.names = FALSE)

##########################################
# (2) E-MTAB (expression atlas - occc and their types) TPM
##########################################

ea_tpm <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/E-MTAB-2770-query-results.tsv", skip = 4, col_names = TRUE)

# NAs fill with 0 - assume no expression 
ea_tpm[ , 3:ncol(ea_tpm)][is.na(ea_tpm[ , 3:ncol(ea_tpm)])] <- 0
# skip first 2 columns as its id only
ea_tpm <- ea_tpm[, -1]
colnames(ea_tpm) <- gsub(" ", "_", colnames(ea_tpm))
colnames(ea_tpm)[colnames(ea_tpm) == "Gene_Name"] <- "Geneid"

sum(duplicated(ea_tpm)) #132
ea_tpm_unique <- ea_tpm %>% group_by(Geneid) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% ungroup()

write.csv(ea_tpm,
          file = "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/E-MTAB_RNAseq_TPM.csv",
          row.names = FALSE)

#######################################
# (3) GSE160692_data- raw to tpm conv
#######################################


GSE160692_data <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE160692_all/GSE160692_OVA_UTE_Raw_Transcripts_GEO.tsv")
table(GSE160692_data$`Gene Symbol`)
duplicate_gene_GSE160692 <- GSE160692_data$`Gene Symbol`[duplicated(GSE160692_data$`Gene Symbol`)]

# only the OVA columns
ova_cols <- grep("^OVA", names(GSE160692_data), value = TRUE)

# group by all non-OVA columns, and average only OVA columns
duplicate_means_GSE160692 <- GSE160692_data %>%
  group_by(across(-all_of(ova_cols))) %>%
  summarise(across(all_of(ova_cols), mean), .groups = "drop")

ova_table <- dplyr::select(duplicate_means_GSE160692, `Gene Symbol`, Start, End, Chromosome, Strand, dplyr::all_of(ova_cols))

# there are gene symbols with dates that start with number migth not match with others
# durign alignment only align the ones match tgt


gene_with_dates <- ova_table %>% filter(grepl("^[0-9]", `Gene Symbol`)) # viewed first

GSE160692_raw_count <- ova_table %>% filter(!grepl("^[0-9]", `Gene Symbol`))

# raw to tpm
gene_id = GSE160692_raw_count$`Gene Symbol`

gr <- GRanges(seqnames = GSE160692_raw_count$Chromosome,
              ranges = IRanges(start = GSE160692_raw_count$Start, end = GSE160692_raw_count$End),
              strand = GSE160692_raw_count$Strand,
              gene_id = GSE160692_raw_count$`Gene Symbol`)


se <- SummarizedExperiment(assays = list(counts = as.matrix(GSE160692_raw_count[, ova_cols])),
                           rowRanges = gr,
                           colData = DataFrame(samples = ova_cols))

counts <- assay(se, "counts")
lengths_kb <- width(rowRanges(se)) / 1000

count2tpm <- function(counts, lengths_kb) {
  rpk <- counts / lengths_kb
  tpm <- t( t(rpk) / colSums(rpk) ) * 1e6
  return(tpm)
}

tpm_GSE160692 <- count2tpm(counts, lengths_kb)
rownames(tpm_GSE160692) <- GSE160692_raw_count$`Gene Symbol`

tpm_GSE160692 <- as.data.frame(tpm_GSE160692)
tpm_GSE160692<- cbind(`Gene Symbol` = rownames(tpm_GSE160692), tpm_GSE160692)
rownames(tpm_GSE160692) <- NULL

colnames(tpm_GSE160692) <- trimws(colnames(tpm_GSE160692))
# tpm_GSE160692 <- tpm_GSE160692 %>% rename(Geneid = `Gene Symbol`) #error on the name
colnames(tpm_GSE160692)[grep("Gene", colnames(tpm_GSE160692))] <- "Geneid"

write.csv(tpm_GSE160692, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE160692_RNA_seq_TPM.csv", row.names = FALSE)


########################################
# (4) tpm file extract ccc GSE189553
########################################


GSE189553_data <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE189553_all/GSE189553_gene_TPM_matrix.txt")

GSE189553_data <- dplyr::select(GSE189553_data, starts_with("gene"), starts_with("CCC"))

GSE189553_data <- dplyr::rename(GSE189553_data, Geneid = gene) 
GSE189553_data<- dplyr::relocate(GSE189553_data, Geneid, .before = 1)

tpm_GSE189553 <- as.data.frame(GSE189553_data)

sum(duplicated(tpm_GSE189553$Geneid)) # 1706

tpm_GSE189553 <- tpm_GSE189553 %>% group_by(Geneid) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

sum(duplicated(tpm_GSE189553$Geneid))

tpm_GSE189553 <- tpm_GSE189553 %>% filter(!grepl("^[0-9]", Geneid))
write.csv(tpm_GSE189553, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE189553_RNAseq_TPM.csv", row.names = FALSE)



###################################
# (5) rpkm to tpm (113 sample)
###################################

rpkm_113_data <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/113_OCCC/RPKM_expdata_113OCCC.xlsx", sheet = 1, col_names = TRUE)

colnames(rpkm_113_data) <- trimws(colnames(rpkm_113_data)) # colname had issue so cleaned
colnames(rpkm_113_data)[1] <- "Geneid"

rpkm_matrix <- as.matrix(rpkm_113_data[, -1])

rpkm2tpm <- function(rpkm_matrix) {
  tpm <- apply(rpkm_matrix, 2, function(x) {
    x / sum(x, na.rm = TRUE) * 1e6})
  return(tpm)}

tpm_matrix <- rpkm2tpm(rpkm_matrix)
tpm_113 <- cbind(Geneid = rpkm_113_data$Geneid, as.data.frame(tpm_matrix))
# no NAs 

sum(duplicated(tpm_113$Geneid))
write.csv(tpm_113, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/tpm113_RNAseq_TPM.csv", row.names = FALSE)

####################################
# (6) rpkm to tpm # paper: Systematic Identification of Characteristic Genes of Ovarian Clear Cell Carcinoma
####################################

rpkm_SD1 <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/ijms-20-04330-s001/Supplmentary Data S1.xlsx", sheet = 1, skip = 2)

rpkm_SD1_OC <- dplyr::select(rpkm_SD1, starts_with("symbol"), starts_with("C"))
rpkm_SD1_Norm <- dplyr::select(rpkm_SD1, starts_with("symbol"), starts_with("N"))

Geneid <- rpkm_SD1_OC[[1]]
rpkm_matrix <- as.matrix(rpkm_SD1_OC[, -1])


rpkm2tpm <- function(rpkm) {
  tpm <- apply(rpkm_matrix, 2, function(x) x / sum(x, na.rm = TRUE) * 1e6)
  return(tpm)}

tpm_matrix <- rpkm2tpm(rpkm_matrix)
tpm_SD1_OC <- cbind(Geneid = Geneid, as.data.frame(tpm_matrix))
tpm_SD1_OC <- as.data.frame(tpm_SD1_OC)

#duplicates
tpm_SD1_OC <- tpm_SD1_OC %>% group_by(Geneid) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

sum(duplicated(tpm_SD1_OC$Geneid))

write.csv(tpm_SD1_OC, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/SD1_RNAseq_TPM.csv", row.names = FALSE)

########################

########################

###################################
# # Kidney
###################################
tpm_ccRCC <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/TCGA-KIRC.star_tpm.tsv")

tpm_ccRCC$Ensembl_ID <- sub("\\..*", "", tpm_ccRCC$Ensembl_ID)


sum(duplicated(tpm_ccRCC$Ensembl_ID))

tpm_ccRCC_mean <- tpm_ccRCC %>% group_by(Ensembl_ID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = tpm_ccRCC_mean$Ensembl_ID,
                   mart = ensembl)

tpm_ccRCC <- tpm_ccRCC_mean %>% left_join(gene_map, by = c("Ensembl_ID" = "ensembl_gene_id"))

# gene symbols have NAs
sum(is.na(tpm_ccRCC)) # [1] 1253
sum(is.na(tpm_ccRCC$hgnc_symbol)) # [1] 1253

# changes to columns
tpm_ccRCC$Geneid <- tpm_ccRCC$hgnc_symbol
tpm_ccRCC$hgnc_symbol <- NULL
tpm_ccRCC$Ensembl_ID <- NULL
tpm_ccRCC <- tpm_ccRCC %>% dplyr:: select(Geneid, everything())

tpm_ccRCC <- as.data.frame(tpm_ccRCC)

#duplicates
tpm_ccRCC <- tpm_ccRCC %>% group_by(Geneid) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

sum(duplicated(tpm_ccRCC$Geneid))

write.csv(tpm_ccRCC, "/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/ccRCC_RNAseq_TPM.csv", row.names = FALSE)

####################################
#proteomics
####################################



####################################

####################################
# ANALYSIS
####################################



##without proteomics##

merge_tpm <- tpm_ccRCC %>% full_join(tpm_SD1_OC, by = "Geneid") %>%
  full_join(tpm_113, by = "Geneid") %>% full_join(tpm_GSE230956, by = "Geneid") %>%
  full_join(tpm_GSE189553, by = "Geneid") %>% full_join(tpm_GSE160692, by = "Geneid") %>% 
  full_join(CCOC_df, by = "Geneid") %>% full_join(ea_tpm_unique, by = "Geneid")

all_sample_names <- colnames(merge_tpm)[-1]


# duplicate check
sapply(list(tpm_ccRCC, tpm_SD1_OC, tpm_113, tpm_GSE230956, tpm_GSE189553, tpm_GSE160692, CCOC_df, ea_tpm_unique), function(df) sum(duplicated(df$Geneid)))
#merge_tpm[is.na(merge_tpm)] <- 0 -> there is non numeric ones
num_cols <- sapply(merge_tpm, is.numeric)
merge_tpm[, -1] <- lapply(merge_tpm[, -1], function(x) { x[is.na(x)] <- 0; x })


# scale
#merge_tpm <- merge_tpm[, sapply(merge_tpm, is.numeric)]
#scaled_tpm <- scale(merge_tpm, center = TRUE, scale = TRUE)

# label and drop id column
kidney_n <- ncol(tpm_ccRCC) - 1
ovary_n <- (ncol(tpm_SD1_OC) - 1) + 
  (ncol(tpm_113) - 1) + 
  (ncol(tpm_GSE230956) - 1) + 
  (ncol(tpm_GSE189553) - 1) + 
  (ncol(tpm_GSE160692) - 1) +
  (ncol(CCOC_df) - 1) +
  (ncol(ea_tpm_unique) - 1)
  
tissue_label <- c(rep("ccRCC", kidney_n), rep("OCCC", ovary_n))

# extract numeric values
exp_matrix <- as.matrix(merge_tpm[, -1])
rownames(exp_matrix) <- merge_tpm$Geneid
#scaled_tpm <- scale(exp_matrix, center = TRUE, scale = TRUE)


#check 0s
table(exp_matrix == 0) 

# 20% non-zero rule
# expressed in at least 20% samples with TPM > 1
expressed_filter <- rowMeans(exp_matrix > 1) >= 0.2  
filtered_exp <- exp_matrix[expressed_filter, ]
# log2 transformation
log_tpm_filtered <- log2(filtered_exp + 1)



#HISTOGRAM

par(mfrow = c(1, 2))
# raw tpm h
hist(as.vector(filtered_exp), breaks = 50, main = "Raw TPM (Filtered Genes)",
     xlab = "TPM",
     col = "skyblue")


# log2(TPM+1) h
hist(as.vector(log_tpm_filtered), breaks = 50, main = "Log2(TPM+1) (Filtered Genes)",
     xlab = "Log2 TPM",
     col = "salmon")


# reseting
par(mfrow = c(1, 1))



# TOP 50 TOGETHER RAW AND LOG

# RAW TPM variance
gene_vars_raw <- apply(filtered_exp, 1, var, na.rm = TRUE)
filtered_gene_vars_raw <- gene_vars_raw[gene_vars_raw > 1 & is.finite(gene_vars_raw)]
n_genes_raw <- min(50, length(filtered_gene_vars_raw))
top_genes_raw <- names(sort(filtered_gene_vars_raw, decreasing = TRUE)[1:n_genes_raw])
top_genes_raw <- intersect(top_genes_raw, rownames(filtered_exp))  # subset
top_raw <- filtered_exp[top_genes_raw, ]


# LOG TPM variance
gene_vars_log <- apply(log_tpm_filtered, 1, var, na.rm = TRUE)
filtered_gene_vars_log <- gene_vars_log[which(gene_vars_log > 1 & is.finite(gene_vars_log))]
n_genes_log <- min(50, length(filtered_gene_vars_log)) # filtered rows the dim not reduced
top_genes_log <- names(sort(filtered_gene_vars_log, decreasing = TRUE)[1:n_genes_log])
top_genes_log <- intersect(top_genes_log, rownames(log_tpm_filtered))  # subset
top_log <- log_tpm_filtered[top_genes_log, ]

#  Z-score normalise each gene (by row) 
top_raw_z <- t(scale(t(top_raw), center = TRUE, scale = TRUE))
top_log_z <- t(scale(t(top_log), center = TRUE, scale = TRUE))


# heatmaps 
ht_raw <- Heatmap(
  cor(t(top_raw_z), method = "spearman", use = "pairwise.complete.obs"),
  name = "(No Log)",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)

ht_log <- Heatmap(
  cor(t(top_log_z), method = "spearman", use = "pairwise.complete.obs"),
  name = "(Log2)",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)

# side by side
draw(ht_raw + ht_log, heatmap_legend_side = "right")

###################
# VIEW SPERATLY RAW
##################


# Heatmap for raw data
ht_raw <- Heatmap(
  cor(t(top_raw_z), method = "spearman", use = "pairwise.complete.obs"),
  name = "ρ (Raw)",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)
draw(ht_raw, heatmap_legend_side = "right")

# with correlation matrix values ########
cor_raw <- cor(t(top_raw_z), method = "spearman", use = "pairwise.complete.obs")
ht_raw <- Heatmap(
  cor(t(top_raw_z), method = "spearman", use = "pairwise.complete.obs"),
  name = "ρ (Raw)",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", cor_raw[i, j]), x, y, gp = gpar(fontsize = 6))
  }
)
draw(ht_raw, heatmap_legend_side = "right")



####################
# VIEW SEPARTLY LOG
###################


# Heatmap for log2-transformed data
ht_log <- Heatmap(
  cor(t(top_log_z), method = "spearman", use = "pairwise.complete.obs"),
  name = "ρ (Log2)",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)
draw(ht_log, heatmap_legend_side = "right")


###### correlation matrix ######### with the values 

cor_mat <- cor(t(top_log_z), method = "spearman", use = "pairwise.complete.obs")

#ha <- HeatmapAnnotation(Histotype = selected_tissue_labels, col = list(Histotype = c("ccRCC" = "blue", "OCCC" = "red")))

ht_log <- Heatmap(
  cor_mat,
  name = " ",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 6))
  }
)

draw(ht_log, heatmap_legend_side = "right")

###################
# common genes
###################
# any common overlaps between genes of log tr and witouth log trasnformed
common_genes <- intersect(top_log_z, top_raw_z)
length(common_genes)  # How many overlap?
print(common_genes) #0

# distance matching instead of exact 
# For each gene in top_genes_log, find closest match in top_genes
# distance metrics (e.g., Levenshtein distance
matches <- sapply(top_log_z, function(g) {
  match_idx <- agrep(g, n_genes_log, max.distance = 0.1, value = TRUE)
  if(length(match_idx) == 0) NA else match_idx
})

print(matches)

##############


#############
#pca
#############

different_papers <- c(
  rep("kidney", ncol(tpm_ccRCC) - 1),
  rep("ovary_SD1_OC", ncol(tpm_SD1_OC) - 1),
  rep("ovary_113", ncol(tpm_113) - 1),
  rep("ovary_GSE230956", ncol(tpm_GSE230956) - 1),
  rep("ovary_GSE189553", ncol(tpm_GSE189553) - 1),
  rep("ovary_GSE160692", ncol(tpm_GSE160692) - 1),
  rep("ovary_CCOC_df", ncol(CCOC_df) - 1),
  rep("ovary_ea_tpm_unique", ncol(ea_tpm_unique) - 1)
)


pca <- prcomp(t(top_log_z), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], col = as.factor(different_papers), pch = 16,
     xlab = "PC1", ylab = "PC2")

legend("bottomright", legend = unique(different_papers), col = 1:length(unique(different_papers)),
       pch = 16, cex = 0.5, pt.cex = 0.8)

summary(pca) # explains the variance
head(pca$rotation) #understand what is the driver genes

# ggplot ver
ggplot(data.frame(pca$x), aes(PC1, PC2, color=different_papers)) +
  geom_point() +
  theme_minimal()


############################

###########
# ----- for log2foldchange only
##########
# filter genes with zero variance in either group
keep_genes <- apply(log_tpm_filtered, 1, function(x) {
  var1 <- var(x[tissue_label == "ccRCC"], na.rm = TRUE)
  var2 <- var(x[tissue_label == "OCCC"], na.rm = TRUE)
  !is.na(var1) && !is.na(var2) && var1 > 0 && var2 > 0
})

log_tpm_filtered2 <- log_tpm_filtered[keep_genes, ]


# mean expression per tissue
mean_Kidney <- rowMeans(log_tpm_filtered2[, tissue_label == "ccRCC"])
mean_Ovary <- rowMeans(log_tpm_filtered2[, tissue_label == "OCCC"])

# log2 fold change: difference of group means 
log2FC <- mean_Kidney - mean_Ovary

# Wilcoxon p-values per gene
pvals_wilcox <- apply(log_tpm_filtered2, 1, function(x) {
  group1_vals <- x[tissue_label == "ccRCC"]
  group2_vals <- x[tissue_label == "OCCC"]
  
  # minimum sample size
  if(length(group1_vals) < 2 || length(group2_vals) < 2) return(NA)
  
  wilcox.test(group1_vals, group2_vals)$p.value
})


# combine 
results_wilcox <- data.frame(
  Gene = rownames(log_tpm_filtered2),
  log2FC = log2FC,
  pvalue = pvals_wilcox
)
# Remove NA p-values - no NAs present
#results_wilcox <- results_wilcox[!is.na(results_wilcox$pvalue), ]

# Adjust p-values (FDR)
results_wilcox$padj <- p.adjust(results_wilcox$pvalue, method = "BH")

cat("Genes with valid Wilcoxon p-values:", nrow(results_wilcox), "\n")
# first volcano plot with many labels
# Enhanced volcano plot with labels 
EnhancedVolcano(results_wilcox,
                lab = results_wilcox$Gene,
                x = 'log2FC',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 3)

# top 10 upregulated in kidney
head(results_wilcox[order(-results_wilcox$log2FC), ], 10)

# top 10 upregulated in ovary
head(results_wilcox[order(results_wilcox$log2FC), ], 10)

###########################

# avoid log(0) issues by replacing zeros with the smallest nonzero padj
min_nonzero <- min(results_wilcox$padj[results_wilcox$padj > 0])
padj_safe <- ifelse(results_wilcox$padj == 0,
                    min_nonzero,
                    results_wilcox$padj)

#  -log10
results_wilcox$neglog10_padj <- -log10(padj_safe)

# histogram of -log10(FDR)
ggplot(results_wilcox, aes(x = neglog10_padj)) +
  geom_histogram(binwidth = 5, color = "white") +
  labs(title = "-log10 Adjusted p-value (FDR) distribution",
       x = "-log10(padj)",
       y = "Count") +
  theme_bw()


# ----- threshold to tweak -----
lfc_cut  <- 1          # |log2FC| >= 1  (~2x)
padj_cut <- 0.05       # FDR cutoff
n_label  <- 20         # how many genes to label
# ------------------------------------

# handle padj==0 
min_nonzero <- min(results_wilcox$padj[results_wilcox$padj > 0], na.rm = TRUE)
padj_safe   <- ifelse(results_wilcox$padj == 0, min_nonzero, results_wilcox$padj)

volc <- results_wilcox %>%
  mutate(
    neglog10_padj = -log10(padj_safe),
    sig = case_when(
      padj <= padj_cut & abs(log2FC) >= lfc_cut ~ "Significant (FDR & |LFC|)",
      TRUE ~ "Not significant"
    )
  )

# labels: top by FDR, then by effect size
top_labs <- volc %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = n_label)

# draw volcano
ggplot(volc, aes(x = log2FC, y = neglog10_padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("Significant (FDR & |LFC|)" = "#d62728",
                                "Not significant" = "grey40")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
  geom_text_repel(data = top_labs,
                  aes(label = Gene),
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.size = 0.2) +
  labs(title = "Volcano plot with top gene labels",
       x = "log2 fold change (ccRCC – OCCC)",
       y = "-log10(FDR)") +
  theme_bw() +
  theme(legend.title = element_blank())
# I USED THIS
# to make it more balanced


# ----- thresholds -----
lfc_cut   <- 1
padj_cut  <- 0.05
n_each    <- 10   # labels per side (total labels = 2*n_each)
# ----------------------

# safe -log10(FDR)
min_nonzero <- min(results_wilcox$padj[results_wilcox$padj > 0], na.rm = TRUE)
padj_safe   <- ifelse(results_wilcox$padj == 0, min_nonzero, results_wilcox$padj)

volc <- results_wilcox %>%
  mutate(
    neglog10_padj = -log10(padj_safe),
    sig = ifelse(padj <= padj_cut & abs(log2FC) >= lfc_cut,
                 "Significant (FDR & |LFC|)", "Not significant")
  )

# pick labels separately on each side:
top_pos <- volc %>%
  filter(log2FC > 0) %>%
  arrange(padj, desc(log2FC)) %>%
  slice_head(n = n_each)

top_neg <- volc %>%
  filter(log2FC < 0) %>%
  arrange(padj, abs(log2FC)) %>%  # smallest padj, large magnitude
  slice_head(n = n_each)

top_labs <- bind_rows(top_pos, top_neg)

# plot
ggplot(volc, aes(x = log2FC, y = neglog10_padj)) +
  geom_point(aes(color = sig), alpha = 0.65, size = 1.2) +
  scale_color_manual(values = c("Significant (FDR & |LFC|)" = "#d62728",
                                "Not significant" = "grey55")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
  geom_text_repel(data = top_labs,
                  aes(label = Gene),
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.size = 0.2) +
  labs(title = "Volcano plot with balanced labels",
       x = "log2 fold change (ccRCC – OCCC)",
       y = "-log10(FDR)") +
  theme_bw() +
  theme(legend.title = element_blank())

###################


###################
##cheks + labels 
stopifnot(ncol(log_tpm_filtered) == length(tissue_label))
tissue_label <- factor(tissue_label, levels = c("ccRCC", "OCCC"))
idx_cc <- which(tissue_label == "ccRCC")
idx_oc <- which(tissue_label == "OCCC")

## overlap of top-variance gene sets
common_genes <- intersect(rownames(top_log_z), rownames(top_raw_z))
message("Common genes: ", length(common_genes))
stopifnot(length(common_genes) > 0)

## build expression matrix for DE (log2(TPM+1))
X <- log_tpm_filtered[common_genes, , drop = FALSE]

## build expression matrix for DE (log2(TPM+1)) 
X <- log_tpm_filtered[common_genes, , drop = FALSE]

## limma DE (OCCC vs ccRCC) 
library(limma)
design <- model.matrix(~ 0 + tissue_label)
colnames(design) <- levels(tissue_label)  # "ccRCC","OCCC"
cont <- makeContrasts(OCCC_vs_ccRCC = OCCC - ccRCC, levels = design)

fit  <- lmFit(X, design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

#  non-colliding name
de <- topTable(fit2, coef = "OCCC_vs_ccRCC", number = Inf, adjust = "fdr")

## thresholds + group means (from TPM scale for interpretability) 
lfc_cut   <- 1
fdr_cut   <- 0.05
tpm_floor <- 1

de$meanTPM_cc <- rowMeans(filtered_exp[rownames(de), idx_cc, drop = FALSE])
de$meanTPM_oc <- rowMeans(filtered_exp[rownames(de), idx_oc, drop = FALSE])

## marker calls 
de$marker_call <- ifelse(
  de$adj.P.Val < fdr_cut & de$logFC >=  lfc_cut & de$meanTPM_oc >= tpm_floor, "OCCC",
  ifelse(de$adj.P.Val < fdr_cut & de$logFC <= -lfc_cut & de$meanTPM_cc >= tpm_floor, "ccRCC", "none")
)

## exports
topN <- 50

markers_OCCC <- subset(de, marker_call == "OCCC")
markers_OCCC <- markers_OCCC[order(markers_OCCC$adj.P.Val, -markers_OCCC$logFC), ]
markers_OCCC_top <- head(markers_OCCC, topN)

markers_ccRCC <- subset(de, marker_call == "ccRCC")
markers_ccRCC <- markers_ccRCC[order(markers_ccRCC$adj.P.Val, markers_ccRCC$logFC), ]
markers_ccRCC_top <- head(markers_ccRCC, topN)

#write.csv(de,                file = "DE_commonGenes_full.csv")
#write.csv(markers_OCCC_top,  file = "Markers_OCCC_top50.csv")
#write.csv(markers_ccRCC_top, file = "Markers_ccRCC_top50.csv")

##heatmap of top panel
panel_genes <- unique(c(rownames(head(markers_OCCC, 25)),
                        rownames(head(markers_ccRCC, 25))))
panel_genes <- intersect(panel_genes, rownames(X))  # safety
panel_mat   <- X[panel_genes, , drop = FALSE]
panel_z     <- t(scale(t(panel_mat)))

library(ComplexHeatmap)
library(circlize)
ha <- HeatmapAnnotation(
  Histotype = tissue_label,
  col = list(Histotype = c(ccRCC = "#3b82f6", OCCC = "#ef4444"))
)
Heatmap(panel_z,
        name = "z",
        top_annotation = ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = TRUE, cluster_columns = TRUE) |>
  draw(heatmap_legend_side = "right")





############################

###########################
# MULTI OMICS
###########################


########proteomics#########

# OCCC - proteomics


pro_L <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/Proteomic_land/source_file_proL.xlsx", sheet = 14, col_names = TRUE)
# 40 samples in CCOC group, Protein expression, also have normal group (n=31)
#sheet 14 - raw protein intensity values from mass spectrometry PulseDIA
pro_L <- pro_L[, colSums(!is.na(pro_L)) > 0]
# metadata
meta <- dplyr::select(pro_L, X, pat_id, Histology8, Group, type)
unique(pro_L$type)
# filter 
pro_L_ccoc   <- pro_L %>% filter(type == "CCOC")
pro_L_normal <- pro_L %>% filter(type == "Normal")

# save 
#saveRDS(prot_ccoc, "proteomics_CCOC_only.rds")
#saveRDS(prot_normal, "proteomics_Normal_only.rds")
nrow(pro_L_ccoc)
nrow(pro_L_normal)

# remove metadata
pro_L_data <- dplyr::select(pro_L_ccoc, -pat_id, -Histology8, -Group, -type)

# pivot long and split uniprotid
pro_L_long <- pro_L_data %>% pivot_longer(-X, names_to = "Protein", values_to = "Intensity") %>%
  separate(Protein, into = c("UniProtID", "Gene"), sep = "_")

# get gene level by mean - avoid duplicates
pro_L_gene <- pro_L_long %>% group_by(Gene, X) %>%
  summarise(Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# pivot back to wide
pro_L_matrix <- pro_L_gene %>%
  pivot_wider(names_from = X, values_from = Intensity) %>%
  column_to_rownames("Gene")
# gene names as a column
pro_L_matrix_out <- pro_L_matrix %>%
  as.data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("Geneid")

write.csv(pro_L_matrix_out, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/Proteomic_land_PROT_TPM.csv", row.names = FALSE)

# log-transform intensities
pro_L_log <- log2(pro_L_matrix + 1)

# z-score normalize per gene
pro_L_z <- t(scale(t(pro_L_log), center = TRUE, scale = TRUE))

# change rownames to colnames and name it Geneid
pro_L_z <- as.data.frame(pro_L_z, check.names = FALSE) %>% mutate(Geneid = rownames(.)) %>%              
  relocate(Geneid, .before = 1) %>% {rownames(.) <- NULL; .}   

####################################
# DEPMAP proteomics

# teh steps I did I did not convert them to TPM thast a naming error - needs change
# what I did: harmonized RPPA → “gene × sample” matrix so still in harmonized RPPA units (often log2-normalized antibody intensities)
##################
# ES2
###################

prot_es2 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/DepMap_Prot/Harmonized_RPPA_CCLE_ES2.csv", col_names = TRUE)
# already normalised by z-score
# get gene symbol
gene_names <- str_extract(colnames(prot_es2)[-1], "(?<=\\().+?(?=\\))")
colnames(prot_es2) <- c("SampleID", gene_names)

# genes in rows
prot_es2_long <- prot_es2 %>% column_to_rownames("SampleID") %>%  t() %>% as.data.frame()
# genes in first col
prot_es2_wide <- prot_es2_long %>% rownames_to_column(var = "Geneid")
write.csv(prot_es2_wide, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/ES2_PROT_TPM.csv", row.names = FALSE)

################
# JHOC-5
################

prot_jhoc5 <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/DepMap_Prot/Harmonized_RPPA_CCLE_JHOC_5.csv", col_names = TRUE)

gene_names_j <- str_extract(colnames(prot_jhoc5)[-1], "(?<=\\().+?(?=\\))")
colnames(prot_jhoc5) <- c("SampleID", gene_names)

# gene in rows
prot_jhoc5_long <- prot_jhoc5 %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame()
# genes in first col
prot_jhoc5_wide <- prot_es2_long %>% rownames_to_column(var = "Geneid")

write.csv(prot_jhoc5_wide, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/JHOC5_PROT_TPM.csv", row.names = FALSE)


#################
# OVMANA
#################

prot_ovmana <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/DepMap_Prot/Harmonized_RPPA_CCLE_OVMANA.csv", col_names = TRUE)

gene_names_o <- str_extract(colnames(prot_ovmana)[-1], "(?<=\\().+?(?=\\))")
colnames(prot_ovmana) <- c("SampleID", gene_names)

# gene in rows
prot_ovmana_long <- prot_ovmana %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame()
# genes in first col
prot_ovmana_wide <- prot_ovmana_long %>% rownames_to_column(var = "Geneid")
write.csv(prot_ovmana_wide, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/OVMANA_PROT_TPM.csv", row.names = FALSE)

################
# combine all proteomics
###############
pro_L_z$Geneid <- as.character(pro_L_z$Geneid)
prot_es2_wide$Geneid <- as.character(prot_es2_wide$Geneid)
prot_jhoc5_wide$Geneid <- as.character(prot_jhoc5_wide$Geneid)
prot_ovmana_wide$Geneid <- as.character(prot_ovmana_wide$Geneid)

all_prot_oc <- pro_L_z %>%
  full_join(prot_es2_wide,   by = "Geneid") %>%
  full_join(prot_jhoc5_wide, by = "Geneid") %>%
  full_join(prot_ovmana_wide, by = "Geneid")


dim(all_prot_oc)
# need to z-score your RNA-seq TPM and MS proteomics values per gene 
# so that all three datasets are on a comparable scale before correlation - no over normalise but standardise basically

all_prot_oc <- all_prot_oc %>% mutate(across(-Geneid, ~ replace(., is.na(.), 0)))



##################
# combine with RNA-seq
##################

# get RNA samples again
#get all rna tpm genes not just 50 so I do z norm again here 
rna_z <- t(scale(t(log_tpm_filtered), center = TRUE, scale = TRUE))
rna_z <- as.data.frame(rna_z, check.names = FALSE)
rna_z[["Geneid"]] <- rownames(rna_z) #row to col
# mv to first col
rna_z <- rna_z[, c("Geneid", setdiff(names(rna_z), "Geneid")), drop = FALSE]
rownames(rna_z) <- NULL
#sum(duplicated(rna_z$Geneid)) #check if 0

# combine
# traceable samples
prot_pref <- all_prot_oc %>% rename_with(~ paste0("PROT_", .), -Geneid)

rna_pref  <- rna_z %>% rename_with(~ paste0("RNA_", .),  -Geneid)

# merge by geneid
# used top_log_z

merged_omics <- full_join(prot_pref, rna_pref, by = "Geneid")
# na : 60268740 replace w 0 - its coz samples do not have those genes
merged_omics <- merged_omics %>% mutate(across(-Geneid, ~ replace(., is.na(.), 0)))


# standardise with z score again
# filter
prot_cols <- grep("^PROT_", names(merged_omics), value = TRUE)
rna_cols  <- grep("^RNA_",  names(merged_omics), value = TRUE)

# with 20% non zero filter
merged_with_filter <- merged_omics %>%
  rowwise() %>%
  mutate(
    prop_prot_nz = mean(c_across(all_of(prot_cols)) != 0),
    prop_rna_nz  = mean(c_across(all_of(rna_cols))  != 0)
  ) %>%
  ungroup() %>%
  filter(prop_prot_nz >= 0.20 & prop_rna_nz >= 0.20) %>%
  dplyr::select(-prop_prot_nz, -prop_rna_nz)

# without non-zero filter
merged_no_filter <- merged_omics

##
z_with_filter <- t(scale(t(merged_with_filter[ , -1]), center = TRUE, scale = TRUE))
z_no_filter   <- t(scale(t(merged_no_filter[ , -1]), center = TRUE, scale = TRUE))


# if skewed apply log transformation
vals <- unlist(z_with_filter[ , -1]) # exclude geneid
vals <- as.numeric(vals)
ggplot(data.frame(value = vals), aes(x = value)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  labs(title = "Distribution of values - 20%", x = "Value", y = "Count")


vals2 <- unlist(z_no_filter[ , -1]) # exclude geneid
vals2 <- as.numeric(vals2)
ggplot(data.frame(value = vals2), aes(x = value)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  labs(title = "Distribution of values NoF", x = "Value", y = "Count")
# no need log trasnformation its already in scale

##################
# corr #numeric matrix
cor_mat_all <- cor(z_with_filter, method = "spearman", use = "pairwise.complete.obs")

# annotations - if false otherwise name as RNA
omics_type <- ifelse(grepl("^PROT_", colnames(cor_mat_all)), "Proteomics", "RNA")
col_anno <- HeatmapAnnotation(Omic = omics_type, 
                              col = list(Omic = c("Proteomics" = "#E69F00", "RNA" = "#56B4E9")))

# HEATMAP
Heatmap(
  cor_mat_all,
  name = "Spearman\nCorrelation",
  top_annotation = col_anno,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE
)



# selecting top 100

# variance per gene 
gene_var <- apply(z_with_filter, 1, var, na.rm = TRUE)
# names
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:100]
# subset mtx
z_top100 <- z_with_filter[top_genes, ]
# spearman correlation 
cor_mat_top100 <- cor(z_top100, method = "spearman", use = "pairwise.complete.obs")
# Annotate 
omics_type <- ifelse(grepl("^PROT_", colnames(cor_mat_top100)), "Proteomics", "RNA")
col_anno <- HeatmapAnnotation(Omic = omics_type,
                              col = list(Omic = c("Proteomics" = "#E69F00", "RNA" = "#56B4E9")))
# heatmap
Heatmap(
  cor_mat_top100,
  name = "Spearman\nCorrelation",
  top_annotation = col_anno,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE
)


################
#
library(enrichR)

#merge bakc the egenid after normalisation
#  from merged_with_filter where col1 is Geneid
stopifnot(exists("merged_with_filter"), "Geneid" %in% colnames(merged_with_filter))

mat <- as.matrix(merged_with_filter[, setdiff(colnames(merged_with_filter), "Geneid"), drop = FALSE])
rownames(mat) <- as.character(merged_with_filter$Geneid)

z_with_filter <- t(scale(t(mat), center = TRUE, scale = TRUE))  # row-wise z
stopifnot(!is.null(rownames(z_with_filter)))


# gene list: e.g. your top 100 variable genes already computed
genes_for_enrichr <- rownames(z_with_filter)[order(apply(z_with_filter,1,var), decreasing=TRUE)][1:100]

dbs_avail <- listEnrichrDbs()

head(dbs)

dbs <- c("Reactome_2022", "KEGG_2021_Human", "GO_Biological_Process_2021")

enr <- enrichr(genes_for_enrichr, dbs)

head(enr$Reactome_2022[, c("Term", "Overlap", "Adjusted.P.value", "Combined.Score")], 20)


library(dplyr); library(tidyr); library(stringr)

rct <- enr$Reactome_2022 %>%
  dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes = Genes) %>%
  separate_wider_delim(Overlap, delim = "/", names = c("k","M"), cols_remove = FALSE) %>%
  mutate(k = as.integer(k), M = as.integer(M), overlap_frac = k / M) %>%
  arrange(Adjusted.P.value)

# Top 15 by adjusted p-value
rct_top <- rct %>% slice_head(n = 15)
rct_top
# Get the hit genes for one pathway (e.g., PI3K/AKT signaling)
rct_df <- as.data.frame(enr$Reactome_2022, stringsAsFactors = FALSE)
ii <- grep("PI3K\\s*/?\\s*AKT", rct_df$Term, ignore.case = TRUE)
if (length(ii) > 0) {
  genes <- unique(unlist(strsplit(rct_df$Genes[ii[1]], "[,;]")))
  genes
}


#################
# DOt plot reactom terms

# tidy tibble
rct <- as.data.frame(enr$Reactome_2022, stringsAsFactors = FALSE) |>
  as_tibble() |>
  dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) |>
  separate_wider_delim(Overlap, "/", names = c("k","M"), cols_remove = FALSE) |>
  mutate(k = as.integer(k), M = as.integer(M),
         neglog10FDR = -log10(Adjusted.P.value)) |>
  arrange(Adjusted.P.value) |>
  slice_head(n = 20) |>
  mutate(Term = factor(Term, levels = rev(Term)))  # order top-to-bottom

ggplot(rct, aes(x = neglog10FDR, y = Term, size = k, fill = Combined.Score)) +
  geom_point(shape = 21, alpha = 0.9) +
  scale_size_area(max_size = 10) +
  labs(x = "-log10(FDR)", y = NULL, size = "Overlap (k)",
       fill = "Combined score", title = "Pathway enrichment in OCCC vs ccRCC") +
  theme_minimal(base_size = 12)

# combiend score: combine significance and effect size
# adj p value: stats sig (plotted as -log10(FDR) - common ranking metirc)

ggplot(rct, aes(x = -log10(Adjusted.P.value), y = Term, size = k, fill = Combined.Score)) +
  geom_point(shape = 21) +
  labs(x = "-log10(FDR)", size = "Overlap (k)", fill = "Combined score")


# otehre than enriched thsi is for down adn up regulated genes
# up genes

# down genes



#keep only variable (non-housekeeping) genes,

#find significant RNA–protein correlations, and

#see whether the relationships are OCCC vs ccRCC–specific

library(broom)
library(matrixStats)
#try try
############
# coampare OCCC and ccRCC
############

# either load houskeepign gene txt or get high variance
gene_var2 <- apply(z_with_filter, 1, var, na.rm = TRUE)
var_thresh <- quantile(gene_var2, 0.50, na.rm = TRUE) # keep top 50% most variable
z_df <- z_with_filter[gene_var2 >= var_thresh, ]

prot_cols <- grep("^PROT_", colnames(z_df), value = TRUE)
rna_cols  <- grep("^RNA_",  colnames(z_df), value = TRUE)

# top 50% variable
gene_var <- apply(z_df[, c(prot_cols, rna_cols)], 1, var, na.rm = TRUE)
var_thresh <- quantile(gene_var, 0.50, na.rm = TRUE)
z_var <- z_df[gene_var >= var_thresh, ]

# tissue label for rna and prot
kidney_l <- ncol(tpm_ccRCC) - 1  # -1 for Geneid column
ovary_l <- (ncol(tpm_SD1_OC) - 1) + 
  (ncol(tpm_113) - 1) + 
  (ncol(tpm_GSE230956) - 1) + 
  (ncol(tpm_GSE189553) - 1) + 
  (ncol(tpm_GSE160692) - 1) +
  (ncol(CCOC_df) - 1) +
  (ncol(ea_tpm_unique) - 1)

rna_labels <- c(rep("ccRCC", kidney_l), rep("OCCC", ovary_l))


prot_ovary_l <- (ncol(pro_L_z) - 1) + 
  (ncol(prot_es2_wide) - 1) + 
  (ncol(prot_ovmana_wide) - 1) +
  (ncol(prot_ovmana_wide) - 1)

prot_labels <- rep("OCCC", prot_ovary_l)


# get names
# RNA sample names, labels
rna_samples <- gsub("^RNA_", "", grep("^RNA_", colnames(z_with_filter), value = TRUE))
meta_rna <- data.frame(SampleID = rna_samples, TumorType = rna_labels)

# proteomics sample names, labels
prot_samples <- gsub("^PROT_", "", grep("^PROT_", colnames(z_with_filter), value = TRUE))
meta_prot <- data.frame(SampleID = prot_samples, TumorType = prot_labels)

# Combine into one metadata frame
meta <- rbind(meta_rna, meta_prot)










#############
#RNA–protein correlation + tumor type func and viz
#############
prot_cols <- grep("^PROT_", colnames(z_with_filter), value = TRUE)
rna_cols  <- grep("^RNA_",  colnames(z_with_filter), value = TRUE)

rna_samples_in_data  <- sub("^RNA_",  "", rna_cols)
prot_samples_in_data <- sub("^PROT_", "", prot_cols)

# OCCC sample IDs
occ_samples <- meta$SampleID[meta$TumorType == "OCCC"]
occ_samples_in_both <- intersect(rna_samples_in_data, prot_samples_in_data)

# Match RNA and PROT columns for only OCCC samples
occ_rna_cols  <- paste0("RNA_",  occ_samples_in_both)
occ_prot_cols <- paste0("PROT_", occ_samples_in_both)


stopifnot(all(occ_rna_cols %in% colnames(z_with_filter)))
stopifnot(all(occ_prot_cols %in% colnames(z_with_filter)))

#  Spearman correlation
spearman_one <- function(x, y) {
  if (length(x) != length(y)) return(tibble(rho = NA, p = NA, n = NA))
  if (sum(!is.na(x) & !is.na(y)) < 3) return(tibble(rho = NA, p = NA, n = NA))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), p = ct$p.value, n = sum(complete.cases(x, y)))
}

results_list <- list()

for (g in merged_with_filter$Geneid) {
  # RNA and PROT values only for OCCC samples
  rna_vals  <- as.numeric(z_with_filter[merged_with_filter$Geneid == g, occ_rna_cols])
  prot_vals <- as.numeric(z_with_filter[merged_with_filter$Geneid == g, occ_prot_cols])
  
  res <- spearman_one(rna_vals, prot_vals) %>% 
    mutate(gene = g, contrast = "OCCC")
  
  results_list[[g]] <- res
}

results_df <- bind_rows(results_list) %>%
  mutate(FDR = p.adjust(p, method = "BH"))

sig_occ <- results_df %>% filter(FDR < 0.05)



setdiff(occ_rna_cols, colnames(z_with_filter))
setdiff(occ_prot_cols, colnames(z_with_filter))


######

##################
# ccRCC proteomics
##################

######
#if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
#if (!require("AnnotationDbi")) BiocManager::install("AnnotationDbi")

prot_renal_raw <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/ccRCC_datasets/TCGA-KIRC.protein.tsv")

prot_renal_raw[is.na(prot_renal_raw)] <- 0

# Strip suffixes like "_pT37T46", "-cleaved", etc.
prot_renal_raw$CleanName <- str_replace(prot_renal_raw$peptide_target, "[-_].*$", "")

# Try mapping protein symbols using ALIAS2EG and then back to SYMBOL
mapped <- mapIds(
  org.Hs.eg.db,
  keys = prot_renal_raw$CleanName,
  column = "SYMBOL",
  keytype = "ALIAS",
  multiVals = "first"
)

# Add mapping to your dataframe
prot_renal_raw$GeneSymbol <- mapped[prot_renal_raw$CleanName] # 487

unmapped <- prot_renal_raw[is.na(prot_renal_raw$GeneSymbol), "peptide_target"]
length(unmapped)
head(unmapped) #161

# Filter to mapped proteins only
prot_renal_mapped <- prot_renal_raw[!is.na(prot_renal_raw$GeneSymbol), ] # 326


colnames(prot_renal_mapped)[colnames(prot_renal_mapped) == "peptide_target"] <- "Geneid"
prot_renal_mapped$Geneid <- toupper(prot_renal_mapped$Geneid) # for capital cases

# Clean the Geneid column by removing anything after first _ or -
prot_renal_mapped$Geneid <- toupper(gsub("[-_].*$", "", prot_renal_mapped$Geneid))

write.csv(prot_renal_mapped, "TCGA-KIRC.protein_with_genes2.csv", row.names = FALSE)






############################
# OTHER ANALYSIS AND PATHWAYS
############################





# label with ccRCC and OCCC - for heatmap with their accompanying genes - top 10

# bar chart of each gene transcript level w/ labels of k vs o - can I ?

# boxplot and heatmap with common and rare next to each other for each gene ( in total 5 gene maybe)

# pathway enriched plot (-logp adj, generatio legend)












############################
# CLINICAL FOR ALL
############################















