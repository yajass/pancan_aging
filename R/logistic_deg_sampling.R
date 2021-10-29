library(wateRmelon)
library(GeneOverlap)
library(pbmcapply)
library(data.table)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(limma)
library(edgeR)
library(TCGAutils)
library(tidyverse)
library(limma)
library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)

# load data
source("code/useful_functions.R")
source("code/02 - survival.R")
survival_data <- survival_data %>%
  filter(Quartiles %in% c(1,4)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))
conf.data <- confounders
conf.data[is.na(conf.data)] = ""
load("../Data/NewTCGAcounts/Full Matrix.RData")
types <- c(aa_cancers, naa_cancers)
# types <- types[types != "KICH"] ## kich has ~ 17 samples/ggroup

## find samples common to rnaseq and clinical
findSamples <- function(cancers = types, counts = RNAseq_counts){
  cancer_age_cont <- survival_data %>% dplyr::filter(type == cancers)
  matched_counts <- matchCountsClinical(filtered_counts = counts, survival_data = cancer_age_cont)
  cancer_counts_cont <- matched_counts[[1]]
  return(colnames(cancer_counts_cont))
}
samps <- pbmclapply(types, findSamples, counts = RNAseq_counts, mc.cores = 10) %>% unlist()

# subset survival data
survival_data <- survival_data %>%
  filter(bcr_patient_barcode %in% samps)

# plot
survival_data %>%
  filter(type %in% total_cancers) %>%
  mutate(class = ifelse(type %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  dplyr::count(class, type, Quartiles) %>%
  mutate(class = factor(class, levels = c("Age-Associated", "Not Age-Associated"))) %>%
  mutate(Quartiles = factor(Quartiles, levels = c("Young", "Old"))) %>%
  ggbarplot(., x = "type", y = "n", fill = "Quartiles", position = position_dodge2(padding = 0.2), stat = "identity", palette = "jco", facet.by = "class", scales = "free") +
  labs(x = NULL, y = "RNA-Seq Sample Size (Primary Tumors)") + theme_pubr(20) + coord_flip() + theme(legend.position = "none")
ggsave("results/figures/fig1/sample_sizes_rnaseq.eps", width = 8, height = 8)


survival_data %>%
  filter(type %in% total_cancers) %>%
  mutate(class = ifelse(type %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  dplyr::count(class, type,) %>%
  mutate(class = factor(class, levels = c("Age-Associated", "Not Age-Associated"))) %>%
  ggplot(., aes(x = class, y = n, fill = class)) +
  geom_bar(stat = "summary", fun = "mean") +
  geom_jitter(size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_ci", width = 0.2, size = 2) +
  labs(x = NULL, y = "RNA-Seq Sample Size (Primary Tumors)") +
  theme_pubr(20) + scale_fill_manual(values = cols) +
  stat_compare_means(method = "wilcox", comparisons = list(c("Age-Associated", "Not Age-Associated")), 
                     label = "p.signif", bracket.size = 1.5, size = 10) +
  theme(legend.position = "none")
ggsave("results/figures/fig1/samp_size_logistic_sample_size.eps", width = 7, height = 7)


