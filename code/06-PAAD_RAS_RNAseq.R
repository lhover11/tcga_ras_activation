# Purpose: perform DE to identify genes differentially expressed between RAS mutated/CN gain and RAS WT

# load libraries
library(tidyverse)
library(edgeR)
library(data.table)
library(limma)
library(ggrepel)

# read in data
meta <- read.csv("../processed/TCGA_PAAD_meta_for_ML.csv")
raw <- read.table("../raw-data/PAAD/TCGA-PAAD.star_counts.tsv", header = TRUE)
data <- read.csv("../processed/TCGA_PAAD_log_norm_counts.csv")
gene_anno <- read.table("../raw-data/PAAD/gencode.v36.annotation.gtf.gene.probemap", skip = 1)

## Look at structure of raw data
raw[1:2, 1:2]

## reduce gene annotation to ensembl id and gene name
gene_anno <- gene_anno[, 1:2]

# Only use samples in both meta data and with rnaseq data
length(names(raw)[!names(raw) %in% meta$sampleid]) # 27 samples in rnaseq and not in meta data
raw <- raw[, c("Ensembl_ID", meta$sampleid)]

# make sure meta data and counts data are in the same order
all(meta$sampleid == names(raw)[2:ncol(raw)])

## go back to raw counts from log2(counts+1)
counts <- as.data.frame(lapply(raw[, 2:ncol(raw)], function(x) 2^x))
counts <- as.data.frame(lapply(counts, function(x) x-1))
genes <- raw$Ensembl_ID


## Find average library size
mean(colSums(counts))/1000000  # 48 million

lib_size <- colSums(counts)
range(lib_size) # 19- 85 mil

## normalize counts
cpm <- cpm(counts)
log <- cpm(counts, log=TRUE)


# Need to remove lowly expressed genes
norm_hist <- log %>%
  as.data.frame() %>%
   pivot_longer(cols = everything(),
                  names_to = "samples",
                  values_to = "norm_counts")
ggplot(norm_hist, aes(x = norm_counts)) +
  geom_histogram(bins = 30, fill = "seagreen", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Log Norm Counts")
## Can see that we have a big peak of genes that aren't expressed in these samples (norm counts ~ -5)

ggplot(norm_hist, aes(x= norm_counts, color = samples)) +
  geom_density() +
  theme_minimal() +
  theme(legend.position = 'none',
        text = element_text(size = 20))

row.names(counts) <- genes
keep <- filterByExpr(counts, group=meta$ras_event)
## function is keeping genes that have 10 read counts in at least 70% of the samples in the non-ras event group (smaller of the 2 groups)
## in this dataset, our mean lib size is 48 mil, 10/48 ~ 0.2.  Function will look for CPM > 0.2 in at least 39 samples
filt_counts <- counts[keep, ] ## 21915 genes


## normalize filtered counts
cpm_filt <- cpm(filt_counts)
log_filt <- cpm(filt_counts, log=TRUE)

norm_hist_filt <- log_filt %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "samples",
               values_to = "norm_counts")
ggplot(norm_hist_filt, aes(x= norm_counts, color = samples)) +
  geom_density() +
  theme_minimal() +
  theme(legend.position = 'none',
        text = element_text(size = 20))


# PCA - still to do


# Create design matrix
group <- meta$ras_event
group <- ifelse(group == 1, "RAS", "WT")
group <- as.factor(group)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

# Set up contrast matrix
contr.matrix <- makeContrasts(
  RASvsWT = RAS-WT, 
  levels = colnames(design))
contr.matrix

# Differential expression analysis
v <- voom(filt_counts, design, plot=TRUE)
v
## variance decreases with increased gene expression

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))
# 4726 genes down in RAS, 3698 genes up in RAS

# get out results
rasvwt <- topTable(efit, coef=1, n=Inf)

## add in gene anno
rasvwt$gene_id <- row.names(rasvwt)
rasvwt <- rasvwt %>% left_join(gene_anno, by = c("gene_id" = "V1"))


#### Plot results ####

# Volcano plot

# Define thresholds for significance
logFC_threshold <- 0.5     # Log fold change threshold, none for now
pvalue_threshold <- 0.05 # p-value threshold

# Create a new column for significance
rasvwt$Significant <- with(rasvwt, ifelse(logFC > logFC_threshold & adj.P.Val < pvalue_threshold, "Up_Significant", "Not Significant"))
rasvwt$Significant <- with(rasvwt, ifelse(logFC < -(logFC_threshold) & adj.P.Val < pvalue_threshold, "Down_Significant", Significant))

# Show top 10 genes with largest positive and negative logFC that are significant
top_positive <- rasvwt %>%
  filter(Significant == "Up_Significant") %>%
  arrange(desc(logFC)) %>%
  slice(1:10)

# Identify top 10 genes with the largest negative logFC
top_negative <- rasvwt %>%
  filter(Significant == "Down_Significant") %>%
  arrange(logFC) %>%
  slice(1:10)

# Combine top positive and negative genes
top_genes <- bind_rows(top_positive, top_negative)


png("../results/plots/PAAD_RASAct_vs_WT_Volcano_Plot.png", width = 800, height = 600)
ggplot(rasvwt, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Not Significant" = "grey", "Up_Significant" = "red", "Down_Significant" = "blue")) +
  theme_minimal() +
  labs(title = "TCGA PAAD Ras Event vs Ras WT",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.text.x= element_text(size= 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank()
        ) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = V2),
                  color = "black",
                  segment.color = "black", segment.size = 0.5,
                  min.segment.length = 0,
                  size = 6, point.padding = 0.3, max.overlaps = 10)
dev.off()
