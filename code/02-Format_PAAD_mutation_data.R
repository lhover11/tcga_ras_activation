# Purpose: prep PAAD TCGA data for machine learning
## 1. read in WES somatic mutation data (SNPs and small INDELs)
## 2. Subset to non-silent mutations
## 3. Identify samples with RAS mutations
## 4. Identify samples with high mutational burden


# load libraries
library(tidyverse)
library(edgeR)
library(data.table)

# setwd
setwd("../code")

# read in data
data <- fread("../raw-data/PAAD/TCGA-PAAD.somaticmutation_wxs.tsv")


#### Remove non-silent mutations ####

# Look at all mutation types
unique(data$effect)

# remove synonymous mutations
data <- data %>%
  filter(!grepl("synonymous_variant", effect))

unique(data$effect)

# other mutation types to remove from data
## intron_variant, non_coding_transcript_exon_variant, 5_prime_UTR_variant, "", downstream_gene_variant, 3_prime_UTR_variant;NMD_transcript_variant,
## mature_miRNA_variant, upstream_gene_variant, 3_prime_UTR_variant, 
## intron_variant;non_coding_transcript_variant, intron_variant;NMD_transcript_variant

remove <- c("intron_variant", "non_coding_transcript_exon_variant", "5_prime_UTR_variant", "", "downstream_gene_variant",
            "3_prime_UTR_variant;NMD_transcript_variant",
            "mature_miRNA_variant", "upstream_gene_variant", "3_prime_UTR_variant",
            "intron_variant;non_coding_transcript_variant", "intron_variant;NMD_transcript_variant")

data <- data %>%
  filter(!effect %in% remove)

unique(data$effect)

#### Identify samples with KRAS, HRAS and NRAS samples ####

ras_muts <- data %>% filter(gene %in% c("KRAS", "NRAS", "HRAS"))
ras_muts <- ras_muts %>% select(sample, gene, Amino_Acid_Change)
ras_muts <- ras_muts %>%
  group_by(sample, gene) %>%
  summarize(aa_change = paste(Amino_Acid_Change, collapse = ","))

# check for duplicated sample ids
any(duplicated(ras_muts$sample))

## what percent of samples have RAS mut?
length(ras_muts$sample)/length(unique(data$sample)) * 100  # 63.5% 
## see if this matches percentages seen in cbioportal:
## KRAS - 63%
## NRAS - 0%
## HRAS - 1%
## Yes, this matches what's seen in cbioportal for PAAD

table(ras_muts$gene) # 107 samples with KRAS mut, 1 sample with HRAS mut

## sanity check for a copule RAS mut samples
data %>% filter(sample == ras_muts$sample[1])
data %>% filter(sample == ras_muts$sample[10])

#### Next label samples with the highest mutations ####

## From the paper: We also removed the samples with the highest mutation burden to remove potential false positives.
## We defined these samples based on five standard deviations above the log10 total non-silent somatic mutation count per sample

# get number of mutations per sample
num_muts <- data %>%
  group_by(sample) %>%
  summarize(count = n())
num_muts$log_nonsilent <- log10(num_muts$count)
std_dev <- sd(num_muts$log_nonsilent)
mean_muts <- mean(num_muts$log_nonsilent)

num_muts$outliers <- ifelse(num_muts$log_nonsilent > (mean_muts + (5*std_dev)), "outlier", "not_outlier")
hist(num_muts$log_nonsilent)

## by author's defintion, only 1 outlier sample with outlier mutation burden, most likely hypermutated sample


#### create meta data for mutation status ####

meta <- data.frame(sample = unique(data$sample))
meta <- meta %>%
  full_join(ras_muts, by = c("sample")) %>%
  full_join(num_muts, by = c("sample"))

names(meta) <- c("sampleid", "ras_mutant", "aa_change", "nonsilent_muts", "log_nonsilent_muts", "mut_outlier")

write.csv(meta, "../processed/TCGA_PAAD_WXS_meta.csv", row.names = FALSE)
