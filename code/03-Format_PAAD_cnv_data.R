# Purpose: prep PAAD TCGA data for machine learning
## 1. read in CN data
## 2. Look at data to make sure it's copy number values and not ratios
## 3. Identify samples with RAS amplifications/copy gains


# load libraries
library(tidyverse)
library(edgeR)
library(data.table)

# read in data
## copy number data
data <- fread("../raw-data/PAAD/TCGA-PAAD.gene-level_ascat3.tsv")


## gene annotation data
anno <- fread("../raw-data/PAAD/gencode.v36.annotation.gtf.gene.probemap")

# Look at data, want to make sure I understand what these values represent.  From documentation, I think they're just copy number at the gene level, not log-transformed or a ratio

## range of values
range(data[, 2:ncol(data)], na.rm = TRUE) # 0 - 94

## put in long format for plotting
long <- data %>%
  pivot_longer(-Ensembl_ID)
table(long$value)

data %>%
  pivot_longer(-Ensembl_ID) %>%
  filter(value <= 10) %>%  # for now filter out the extremely high values
  ggplot(aes(value)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(0, 10, by = 1))
# vast majority of values are 2, this data is mostly likely copy number values (2 = diploid)  

# plot log-transformed version since we have some high outlier values
data %>%
  pivot_longer(-Ensembl_ID) %>%
  ggplot(aes(log2(value))) +
  geom_histogram()


# Look at RAS genes (HRAS, NRAS, KRAS)

## add in gene annotation info
anno <- anno %>% select(id, gene)
data <- anno %>% left_join(data, by = c("id" = "Ensembl_ID"))

ras <- data %>%
  filter(gene %in% c("KRAS", "HRAS", "NRAS"))

ras <- ras %>%
  pivot_longer(-c(id, gene))

hist(ras$value)
# most values are low (<=2), but we have some samples amplified RAS genes

# label samples with a copy number > 4 as RAS amplified
ras$amp <- ifelse(ras$value >= 4, "amp", "non")
table(ras$amp) # 68 amplified samples, that seems like a lot
# maybe some of those are duplicate samples, have more than 1 ras gene amplified
length(unique(ras$name[ras$amp == "amp"])) # 40 samples with amplified ras gene, still seems like a lot

# differentiate between samples with CN = 4 and those > 4
ras <- ras %>%
  mutate(amp = case_when(
    value < 4 ~ "no_gain",
    value == 4 ~ "dup",
    value >4 ~ "amp"
  ))
table(ras$amp)  # 43 with cn = 4, 25 with cn > 4, start analysis with cn > 4

# sanity check, make sure the KRAS amplified samples in cbioportal are captured in our amplified samples
ras %>% filter(grepl("TCGA-HZ-8636|TCGA-IB-7886|TCGA-IB-7890|TCGA-HZ-7289", name)) %>%
  arrange(name)
## yes, all samples in cbioportal with KRAS amp show KRAS amplification in our data
## TCGA-HZ-8636 = KRAS CN = 5
## TCGA-IB-7886 = KRAS CN = 7
## TCGA-HZ-7289 = KRAS CN = 21
## TCGA-IB-7890 = KRAS CN = 9 

# Write out meta data for these samples to show amp status
meta <- ras %>% select(gene, name, amp)
meta <- meta %>% filter(amp != "no_gain") # filter to copy gains
meta <- meta %>%
  group_by(name) %>%
  summarize(gene = paste(gene, collapse = ","),
            amp = paste(amp, collapse = ","))

# add in all sample ids so we can see all TCGA IDs which we have CNV data for
all_samples <- data.frame(name = unique(long$name))
meta <- full_join(meta, all_samples, by = c("name"))

names(meta) <- c("sampleid", "ras_gene", "cn_event")

write.csv(meta, "../processed/TCGA_PAAD_RAS_Amp.csv", row.names = FALSE)
