library(tidyverse)
library(pheatmap)
library(pbmcapply)
library(TCGAbiolinks)
library(TCGAutils)
library(GSVA)
library(tidymodels)
library(ggsci)
library(fgsea)

source("code/useful_functions.R")
source("code/02 - survival.R")
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/NewTCGAcounts/Full Matrix.RData")
senesc_pathways <- gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/senesccence_gene_sets.gmt")

# get gene lengths and prep data
gene_info <- TCGAbiolinks::geneInfoHT
common <- intersect(rownames(RNAseq_counts), rownames(gene_info)) # 23186
RNAseq_counts <- RNAseq_counts[rownames(RNAseq_counts) %in% common, ]
RNAseq_counts <- RNAseq_counts[order(rownames(RNAseq_counts)), ]
gene_info <- gene_info[rownames(gene_info) %in% common, ]
gene_info <- gene_info[order(rownames(gene_info)), ]
all(rownames(RNAseq_counts) == rownames(gene_info))

# convert to TPM - https://support.bioconductor.org/p/91218/
RNAseq_counts <- RNAseq_counts / gene_info$geneLength
tpm.mat <-  data.frame(t( t(RNAseq_counts) * 1e6 / colSums(RNAseq_counts) ))
colnames(tpm.mat) <- gsub("\\.","-",colnames(tpm.mat))
rm(RNAseq_counts)

# subset tumor types of interest
survival_data <- survival_data %>%
  filter(type %in% c(aa_cancers, naa_cancers))
match_datasets <- matchCountsClinical(filtered_counts = tpm.mat, survival_data = survival_data, onlyPrimary = TRUE, primaryAndNormal = FALSE)
survival_data <- match_datasets$survival_data
tpm.mat <- match_datasets$filtered_counts
all(colnames(tpm.mat) == survival_data$bcr_patient_barcode); rm(match_datasets)

# convert to gene symbol
tpm.mat <- getGeneNamesFromENSELBL(tpm.mat)
saveRDS(tpm.mat, "results/rds/tpm_matrix.rds")
tpm.mat <- readRDS("results/rds/tpm_matrix.rds")

# run ssgsea
set.seed(20)
enrich.mat <- GSVA::gsva(expr = as.matrix(tpm.mat), gset.idx.list = senesc_pathways, kcdf = "Gaussian",
                         method = "ssgsea", verbose = TRUE, parallel.sz = 14) # does not work - BPPARAM = MulticoreParam(workers = 10, progressbar = verbose, RNGseed = 2019)
saveRDS(enrich.mat, "results/rds/enrichment_scores_senescence.rds"); rm(tpm.mat)
enrich.mat <- data.frame(readRDS("results/rds/enrichment_scores_senescence.rds"))
melt.mat <- enrich.mat %>% 
  rownames_to_column("pathway") %>%
  reshape2::melt("pathway") %>%
  mutate(variable = gsub("\\.","-",variable)) %>%
  dplyr::rename(bcr_patient_barcode = variable) %>%
  inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis, type, Quartiles)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", NA)),
         pathway = gsub("\\_", " ", pathway)) %>%
  drop_na()

t.res <- melt.mat %>%
  group_by(type, pathway)%>%
  rstatix::t_test(value ~ Quartiles) %>%
  ungroup() %>% group_by(type) %>% 
  rstatix::adjust_pvalue("p", "fdr") ## %>%
  # filter(fdr < 0.05)

library(pals)
set.seed(19)
cols = setNames(sample(kelly(), 16),
                t.res$type %>% unique)
  
t.res %>%
  mutate(direction = ifelse(statistic < 0, "Young", "Old"),
         class = ifelse(type %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  ggplot(., aes(x = pathway, y = -log10(fdr), size = -log10(fdr), shape = direction, color = type)) +
  facet_wrap(~class, nrow = 2) + expand_limits(y = 0) +
  geom_point(stroke = 2) + geom_hline(yintercept = -log10(0.05), lty = 2, size = 2) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) +
  scale_shape_manual(values = c(1, 20))  +
  scale_size_continuous(range = c(0.1,10)) +
  scale_color_manual(values = cols) +
  labs(x = NULL, y = "-log FDR", shape = "Enrichment", size = "-log FDR", color = NULL) +
  theme_pubclean(30) +
  guides(color = guide_legend(override.aes = list(size=8))) +
  guides(shape = guide_legend(override.aes = list(size=8)))
ggsave("results/figures/fig3/senescence_aging_ssgsea.eps", width = 35, height = 15)
