library(pheatmap)
library(gtools)
library(reshape2)
library(openxlsx)
library(TCGAutils)
library(matrixStats)
library(tidyverse)

source("code/useful_functions.R")
source("code/02 - survival.R")
# aa_maf <- readRDS("results/rds/aa_maf.rds")
# naa_maf <- readRDS("results/rds/naa_maf.rds")

# Subset age data
survival_data <- survival_data %>%
  filter(Quartiles %in% c(1,4),
         type %in% aa_cancers)
survival_data <- data.frame(TCGA.Participant.Barcode = survival_data$bcr_patient_barcode,
                            Quartile = survival_data$Quartiles,
                            Cancer = survival_data$type)

# Load ddr data
# Genomic and Molecular Landscape of DNA Damage Repair Deficiency across The Cancer Genome Atlas
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5961503/
gene_mutations <- read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA_DDR_Data_Resources/GeneMutations.tsv")
genes <- t(read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA_DDR_Data_Resources/Genes.tsv"))[1,]
sampleNames <- t(read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA_DDR_Data_Resources/Samples.tsv"))[1,]

# Set names for ddr data
rownames(gene_mutations) <- genes
colnames(gene_mutations) <- sampleNames
colnames(gene_mutations) <- substring(colnames(gene_mutations), 1, 12)
gene_mutations <- as.data.frame(t(gene_mutations))
gene_mutations <- rownames_to_column(gene_mutations, "TCGA.Participant.Barcode")

# Merge ddr with age data and order by cancer type
gene_mutations <- merge(survival_data, gene_mutations, "TCGA.Participant.Barcode")
gene_mutations <- gene_mutations[order(gene_mutations$Cancer), ]

# Prep for heatmap
gene_mutations <- na.omit(gene_mutations)
gene_mutations$Quartile <- ifelse(gene_mutations$Quartile == 1, "Young", "Old")
rownames(gene_mutations) <- NULL
gene_mutations <- column_to_rownames(gene_mutations, "TCGA.Participant.Barcode")
gene_mutations <- gene_mutations[with(gene_mutations, order(Cancer, Quartile)), ]    # sort
forHeatmap <- gene_mutations[,3:ncol(gene_mutations)]
forAnnotations <- gene_mutations[,c(1,2)]

# Visualize
pheatmap(forHeatmap[gene_mutations$Cancer %in% c( "UCEC"),], scale = "none", cluster_rows = F,
         annotation_row = forAnnotations[gene_mutations$Cancer %in% c("UCEC"),],
         color = c("white", "black"),show_rownames = F, show_colnames = F, legend = F,
         annotation_colors = list(Quartile = c(Old = pal_jco()(2)[2], Young = pal_jco()(2)[1])),
         filename = "results/figures/fig4/ucec_ddr_binary_heatmap.jpeg")
dev.off()
