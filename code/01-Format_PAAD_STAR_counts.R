# Purpose: prep PAAD TCGA data for machine learning
## 1. read in STAR counts downloaded from xenabrowser
## 2. anti-log values to get back to raw counts
## 3. CPM normalize
## 4. log2 normalized counts & write out for DE
## 5. subset top 8000 most variable genes for ML & write out


# load libraries
library(tidyverse)
library(edgeR)
library(data.table)

# setwd
setwd("../code")

# number of gene annotation columns
anno_cols <- 6

# read in data
data <- fread("../raw-data/PAAD/TCGA-PAAD.star_counts.tsv")

# gene mapping
genes_anno <- fread("../raw-data/PAAD/gencode.v36.annotation.gtf.gene.probemap")

# meta data
meta <- fread("../raw-data/PAAD/TCGA-PAAD.clinical.tsv")


# check we have the expected 183 samples
dim(data) # 60660 rows, 184 columns -- matches what we would expect


# check for any normal samples in the RNAseq data
norms <- meta$sample[meta$sample_type.samples == "Solid Tissue Normal"]
names(data)[names(data) %in% norms]

## remove normal tissue samples from data
data <- as.data.frame(data)
data <- data[, !names(data) %in% norms]
dim(data)

# anti-log counts
## methods from xena: log2(count+1)
genes <- data$Ensembl_ID
counts <- as.data.frame(lapply(data[, 2:ncol(data)], function(x) 2^x))
## remove pseudocount
counts <- as.data.frame(lapply(counts, function(x) x-1))

# Reshape to long format for plotting
counts_hist <- pivot_longer(counts,
                            cols = everything(),
                            names_to = "samples",
                            values_to = "raw_counts")

# Plot histogram of raw counts
ggplot(counts_hist, aes(x = raw_counts)) +
  geom_histogram(bins = 30, fill = "seagreen", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Raw Counts")

# calc cpm
log <- as.data.frame(cpm(counts, log = TRUE, prior.count = 1))


norm_hist <- pivot_longer(log,
                            cols = everything(),
                            names_to = "samples",
                            values_to = "norm_counts")

# Plot histogram of norm counts
ggplot(norm_hist, aes(x = norm_counts)) +
  geom_histogram(bins = 30, fill = "seagreen", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Log Norm Counts")

## obviously have many genes without expression in the PAAD dataset


# from the paper: subset the gene expression matrix to the 8,000 most variably expressed genes by median absolute deviation (MAD)
# Calc MAD without scaling

mad <- apply(log, 1, function(x) mad(x, constant = 1))
range(mad)

# add genes back in 
log$genes <- genes
log <- left_join(log, genes_anno, by = c("genes" = "id"))
log <- log %>% select(genes, gene, chrom, chromStart, chromEnd, strand, everything())

# add in MAD and keep top 8000 most variable
log$mad <- mad
top_8000 <- log %>% top_n(8000, wt = mad)

## look at top most variable genes
top_8000 %>% arrange(-mad) %>% head()

## check for KRAS, NRAS, HRAS
top_8000 %>% filter(grepl("KRAS|HRAS|NRAS", gene)) # not in top most variable genes
log %>% filter(grepl("KRAS|HRAS|NRAS", gene)) # sanity check to make sure genes are there


# Plot histogram of 8000 most variable
top_8000 %>% pivot_longer(cols = contains('TCGA'),
               names_to = "samples",
               values_to = "norm_counts") %>%
  ggplot(aes(x = norm_counts)) +
  geom_histogram(bins = 30, fill = "seagreen", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Most Variable Genes")


# Write out top most variable genes
write.csv(top_8000, "../processed/TCGA_PAAD_top8000_variable.csv", row.names = FALSE)

# Write out log normalized data
write.csv(log, "../processed/TCGA_PAAD_log_norm_counts.csv", row.names = FALSE)
