# updated dataset and analysis script 

library(tidyverse)
library(readxl)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)

#if (!require("BiocManager", quietly = TRUE))
 #install.packages("BiocManager")
#BiocManager::install("GenomicRanges", force = TRUE)
#BiocManager::install("SummarizedExperiment", force = TRUE)

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



###########################
# PROTEOMICS SECTION
###########################


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
#pro_L_log <- log2(pro_L_matrix + 1)

# z-score normalize per gene
#pro_L_z <- t(scale(t(pro_L_log), center = TRUE, scale = TRUE))

# change rownames to colnames and name it Geneid
#pro_L_z <- as.data.frame(pro_L_z, check.names = FALSE) %>% mutate(Geneid = rownames(.)) %>%              
  #relocate(Geneid, .before = 1) %>% {rownames(.) <- NULL; .}   

####################################
# DEPMAP proteomics

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















