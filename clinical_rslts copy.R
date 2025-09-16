# clinical heterogeneity analysis
library(tidyverse)

clinical_56 <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE230956_all/Clinical_56.xlsx", col_names = TRUE)
colnames(clinical_56) <- gsub(" ", "_", colnames(clinical_56))

ovary_data <- clinical_56[clinical_56$Tissue_type == "Ovary", ]

