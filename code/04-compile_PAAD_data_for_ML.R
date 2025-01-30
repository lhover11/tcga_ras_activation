# Purpose: prep PAAD TCGA data for machine learning
## Compile data from mutations and CN to identify the RAS mutants for ML


# load libraries
library(tidyverse)
library(edgeR)
library(data.table)


# read in data
## rnaseq (top 8000 most variable genes)
rna <- fread("../processed/TCGA_PAAD_top8000_variable.csv")

## muts
mut <- fread("../processed/TCGA_PAAD_WXS_meta.csv")

## cnv
cn <- fread("../processed/TCGA_PAAD_RAS_Amp.csv")


# combine muts and cn data into meta data

## convert cn to a simple yes/no status for RAS amp events (cn > 4)
cn$ras_amp <- 0
cn$ras_amp <- ifelse(grepl("amp", cn$cn_event), 1, 0)

cn <- cn %>% select(sampleid, ras_amp)

## convert muts to simple yes/no status for RAS mutation
mut$ras_mut <- 0
mut$ras_mut <- ifelse(!is.na(mut$ras_mutant), 1, 0)

## sample to remove due to hypermutated status
remove <- mut$sampleid[mut$mut_outlier == "outlier"]
mut <- mut %>% filter(sampleid != remove)

mut <- mut %>% select(sampleid, ras_mut)


## combine mut and cn status, only keep samples in both the mut and copy number datasets
combine <- cn %>% inner_join(mut, by = c("sampleid"))

## only keep samples with rnaseq data
combine$sampleid <- gsub("-", ".", combine$sampleid)
combine <- combine %>% filter(sampleid %in% names(rna))

## create ras event column for if sample has ras mut or ras amp
combine$ras_event <- ifelse(combine$ras_amp == 1 | combine$ras_mut == 1, 1, 0)
### sanity check
valid_condition <- combine$ras_event == 1 & (combine$ras_amp == 1 | combine$ras_mut == 1)
all(valid_condition | combine$ras_event == 0) # TRUE

combine <- combine %>% select(sampleid, ras_event)


# Prep RNAseq data for ML, transform so rows are samples and columns are genes
## subset df
rna <- rna %>% select(genes, contains("TCGA"))  # use ensembl id to avoid duplicated gene names

## transform df
rna_ml <- rna %>%
  pivot_longer(
    cols = starts_with("TCGA"),
    names_to = "sampleid",
    values_to = "value"
  )

rna_ml <- rna_ml %>%
  pivot_wider(
    names_from = genes,
    values_from = value
  )
dim(rna_ml)

## subset to samples with cn and mut status
rna_ml <- rna_ml %>% filter(sampleid %in% combine$sampleid)

# Write out files for machine learning step
write.csv(combine, "../processed/TCGA_PAAD_meta_for_ML.csv", row.names = FALSE)
write.csv(rna_ml, "../processed/TCGA_PAAD_top8000_var_for_ML.csv", row.names = FALSE)