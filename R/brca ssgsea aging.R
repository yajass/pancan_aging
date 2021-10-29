###########################################################################################################
###########################################################################################################
####                                        Load functions                                             ####
####                                           and data                                                ####
###########################################################################################################
###########################################################################################################


# Converts ENSEMBL IDs (row names) to unique gene symbols
getGeneNamesFromENSELBL <- function(countData){
  library(org.Hs.eg.db)
  # rownames(countData) <- gsub("\\..*","",rownames(countData))
  countData[,ncol(countData)+1] <- mapIds(org.Hs.eg.db, keys = gsub("\\..*","",rownames(countData)),
                                          keytype = "ENSEMBL", column="SYMBOL", "first")
  countData <- na.omit(countData)
  countData <- countData[-which(duplicated(countData[,ncol(countData)])),]
  rownames(countData) <- countData[,ncol(countData)]
  countData <- countData[,-ncol(countData)]
  return(countData)
}

# Download data
###########################################################################################################
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  legacy = FALSE,
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Solid Tissue Normal")
GDCdownload(query)
temp <- GDCprepare(query)
normal_data <- data.frame(assay(temp))
normal_data <- getGeneNamesFromENSELBL(normal_data)
colnames(normal_data) <- gsub('\\.', '-', substr(colnames(normal_data), 1, 12))

# Survival data
###########################################################################################################
source("code/useful_functions.R")
source("code/02 - survival.R")
brca_samples <- survival_data$bcr_patient_barcode[survival_data$type == "BRCA"]

# Dont use GTEx sig. Use mSigDB signature instead
# # Load signature
# ###########################################################################################################
# res_df <- readRDS("/Users/Yajas/Documents/Elemento/AgeTCGA-master/res_df.rds")
# aging_up <- lapply(res_df, function (u) rownames(u)[u$logFC > 1 & u$adj.P.Val < 0.05])
# aging_down <- lapply(res_df, function (u) rownames(u)[u$logFC < -1 & u$adj.P.Val < 0.05])
# breast_sig <- list(aging_up = aging_up$Breast, aging_down = aging_down$Breast)


# aging signature obtained from
# http://software.broadinstitute.org/gsea/msigdb/cards/LY_AGING_OLD_UP
# Gene Set: LY_AGING_OLD_UP
# Use this signature instead of GTEx
agingdown <- c("AC027237.1","ATR","BARD1","CCNA2","CCNB1","CCNF","CDC20","CDC25B","CDH11","CDK4","CENPA","CENPF","CKAP5","CKS1B","CSE1L","CTSC","CXCL8","DDX39A","FBL","FBN2","FOXM1","H2AFX","H2AFZ","HAS2","HMGB2","HMGN2","HSD17B10","KIF11","KIF14","KIF2C","MCM2","MYBL2","NAE1","NASP","NUP88","PAFAH1B1","PARP1","PCNA","PKMYT1","PLK1","POSTN","PPP1CC","PSMA2","PSMA3","PSMC2","PSMC6","PSMD11","PSMD12","PTGS2","RANBP1","SAFB","SERPINB2","TGFBR2","TYMS","UBE2C","UGCG")
agingup <- c("COMP","CRYAB","CST6","FMOD","HTRA1","MMP12","PTGS1")
senesc_pathways <- fgsea::gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/senesccence_gene_sets.gmt")
breast_sig <- senesc_pathways


###########################################################################################################
###########################################################################################################
####                                        Normal ssGSEA                                              ####
###########################################################################################################
###########################################################################################################


# Select data
###########################################################################################################
survival_data <- survival_data[survival_data$type %in% "BRCA", ] %>%
  filter(Quartiles %in% c(1, 4))

normal_data <- normal_data[, colnames(normal_data) %in% survival_data$bcr_patient_barcode]

# Run ssgsea
###########################################################################################################
library(GSVA)
set.seed(2019)
normal_ssgsea <- data.frame(t(gsva(as.matrix(normal_data), breast_sig, "Poisson", "ssgsea")))

# Combine dataframes
###########################################################################################################
normal_ssgsea <- rownames_to_column(normal_ssgsea, "bcr_patient_barcode")
normal_ssgsea <- merge(normal_ssgsea, survival_data, "bcr_patient_barcode")
normal_ssgsea$Dataset <- "TCGA-BRCA Normal"
normal_ssgsea <- data.frame(normal_ssgsea[, c(1:13,48, 51)])
normal_ssgsea$Quartiles <- ifelse(normal_ssgsea$Quartiles == 1, "Young", "Old")


###########################################################################################################
###########################################################################################################
####                                        Cancer ssGSEA                                              ####
###########################################################################################################
###########################################################################################################

# Load data
###########################################################################################################
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/NewTCGAcounts/Full Matrix.RData")
filtered_counts <- RNAseq_counts
filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Solid Tumor" |
                                        TCGAbiospec(colnames(filtered_counts))$sample_definition == "Primary Blood Derived Cancer - Peripheral Blood")]
biospec_table <- TCGAbiospec(colnames(filtered_counts))
biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
filtered_counts <- filtered_counts[,colnames(filtered_counts)[order(biospec_table$plate, decreasing = TRUE)]]

if (sum(duplicated(biospec_table$submitter_id)) > 0){
  filtered_counts <- filtered_counts[, !duplicated(biospec_table$submitter_id)]
}

# change column names to match metadata
colnames(filtered_counts) <- substr(colnames(filtered_counts), 1, 12)


# Select data
###########################################################################################################
filtered_counts <- filtered_counts[,colnames(filtered_counts) %in%  brca_samples]
filtered_counts <- getGeneNamesFromENSELBL(filtered_counts)

# Run ssgsea
###########################################################################################################
set.seed(2019)
cancer_ssgsea <- data.frame(t(gsva(as.matrix(filtered_counts), breast_sig, "Poisson", "ssgsea")))

# Combine dataframes
###########################################################################################################
cancer_ssgsea <- rownames_to_column(cancer_ssgsea, "bcr_patient_barcode")
cancer_ssgsea <- merge(cancer_ssgsea, survival_data, "bcr_patient_barcode")
cancer_ssgsea$Dataset <- "TCGA-BRCA Cancer"
cancer_ssgsea <- data.frame(cancer_ssgsea[, c(1:14,48, 51)])
cancer_ssgsea$Quartiles <- ifelse(cancer_ssgsea$Quartiles == 1, "Young", "Old")


###########################################################################################################
###########################################################################################################
####                                        METABRIC ssGSEA                                            ####
###########################################################################################################
###########################################################################################################

# Load TCGA data
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/R_scripts/R/Fig.1/clinical_attributes_os_stage.R")
survival_data <- survival_data[survival_data$type == "BRCA", ]
survival_data <- survival_data[survival_data$Quartiles %in% c(1, 4), ]
tcga_q1_max <- max(survival_data$age_at_initial_pathologic_diagnosis[survival_data$Quartiles == 1])
tcga_q4_min <- min(survival_data$age_at_initial_pathologic_diagnosis[survival_data$Quartiles == 4])

# Convert metabric into young and old
metabric_clinical <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_clinical_patient.txt",
                              sep = "\t", header = T, stringsAsFactors = F)
metabric_clinical <- metabric_clinical[-c(1,2), ]
metabric_clinical$Quartile = ifelse(metabric_clinical$Age.at.Diagnosis <= tcga_q1_max, "Young", 
                                    ifelse(metabric_clinical$Age.at.Diagnosis >= tcga_q4_min, "Old", NA))
metabric_clinical <- metabric_clinical[!is.na(metabric_clinical$Quartile), ]

# Load metabric gene expression matrix
metabric_expression <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_expression.txt", 
                                sep = "\t", header = T, stringsAsFactors = F)
metabric_expression <- metabric_expression[, -2]
metabric_expression <- column_to_rownames(metabric_expression, "Hugo_Symbol")
colnames(metabric_expression) <- gsub('\\.', '-', colnames(metabric_expression))

# identify common cases
common_cases <- intersect(colnames(metabric_expression), metabric_clinical$X.Patient.Identifier)

# Keep common samples
metabric_expression <- metabric_expression[,common_cases]
metabric_clinical <- metabric_clinical[metabric_clinical$X.Patient.Identifier %in% common_cases, ]

# ssgsea
set.seed(2019)
metabric_gsea <- data.frame(t(gsva(as.matrix(metabric_expression), breast_sig, "Poisson", "ssgsea")))
metabric_gsea <- rownames_to_column(metabric_gsea, "X.Patient.Identifier")
metabric_gsea$Dataset <- "METABRIC"

# create a clinical subset and combine
subset_clinical <- metabric_clinical[, c("X.Patient.Identifier", "Quartile")]
metabric_gsea <- merge(metabric_gsea, subset_clinical)
colnames(metabric_gsea)[1] <- "bcr_patient_barcode"
colnames(metabric_gsea)[16] <- "Quartiles"
# metabric_gsea <- metabric_gsea[, c("bcr_patient_barcode",
#                                    "aging_up",
#                                    "aging_down",
#                                    "Quartiles",
#                                    "Dataset")]
metabric_gsea$Quartile <- factor(metabric_gsea$Quartiles, levels = c("Young", "Old"))



###########################################################################################################
###########################################################################################################
####                                        GTEx Breast ssGSEA                                         ####
###########################################################################################################
###########################################################################################################

# Load GTEx data
###########################################################################################################
meta1 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
meta2 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t")
gtex_data <- read.table("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/Tissues/GTEx_Breast - Mammary Tissue.txt")

# Subset the data as per young and old groups
###########################################################################################################
# match column names with meta2
colnames(gtex_data) <- gsub('\\.', '-', colnames(gtex_data))
colnames(gtex_data) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", colnames(gtex_data))

# convert transcript rownames to ensembl
gtex_data <- gtex_data[!duplicated(gsub("\\..*","",rownames(gtex_data))), ]
rownames(gtex_data) <- gsub("\\..*","",rownames(gtex_data))

age_meta <- meta2[meta2$SUBJID %in% colnames(gtex_data), ]
age_meta$AGE <- ifelse(age_meta$AGE == "20-29", "Young",
                       ifelse(age_meta$AGE == "30-39", "Young",
                              ifelse(age_meta$AGE == "40-49", "Young",
                                     ifelse(age_meta$AGE == "50-59", NA, "Old"))))
age_meta <- na.omit(data.frame(SUBJID = age_meta$SUBJID, AGE = age_meta$AGE))
gtex_data <- gtex_data[, colnames(gtex_data) %in% age_meta$SUBJID]
gtex_data <- getGeneNamesFromENSELBL(gtex_data)


# Run ssGSSEA and clean up results
###########################################################################################################
set.seed(2019)
gtex_ssgsea <- data.frame(t(gsva(as.matrix(gtex_data), breast_sig, "Poisson", "ssgsea")))
gtex_ssgsea <- rownames_to_column(gtex_ssgsea, "SUBJID")
gtex_ssgsea$Dataset <- "GTEx Breast"
gtex_ssgsea <- merge(gtex_ssgsea, age_meta, "SUBJID")

colnames(gtex_ssgsea)[1] <- "bcr_patient_barcode"
colnames(gtex_ssgsea)[16] <- "Quartiles"
# gtex_ssgsea <- gtex_ssgsea[, c("bcr_patient_barcode",
#                                "aging_up",
#                                "aging_down",
#                                "Quartiles",
#                                "Dataset")]




###########################################################################################################
###########################################################################################################
####                                        Merge Normal and Cancer                                    ####
####                                               Results                                             ####
###########################################################################################################
###########################################################################################################

dat <- bind_rows(metabric_gsea, gtex_ssgsea) %>% 
  bind_rows(normal_ssgsea) %>%
  bind_rows(cancer_ssgsea)
dat$Type = ifelse(dat$Dataset %in% c('GTEx Breast', 'TCGA-BRCA Normal'), 'Normal', 'Cancer')
dat$Quartiles <- factor(dat$Quartiles, levels = c('Young', 'Old'))
colnames(dat)

st = dat %>% group_by(Dataset) %>% rstatix::wilcox_test(GO_AGING ~ Quartiles) %>% rstatix::adjust_pvalue(method = 'fdr') %>% rstatix::add_significance(p.col = 'p.adj')
st
p <- ggplot(dat, aes(x = Quartiles, y = GO_AGING, fill = Quartiles)) + geom_violin() + geom_boxplot(width = 0.2) + 
  facet_wrap(~Type + Dataset, ncol = 2, scales = 'free') + 
  scale_fill_jco() + labs(x = NULL, y = 'ssGSEA Score') + theme_pubr(25, legend = 'none') +
  stat_compare_means(method = 'wilcox', label = 'p.signif', label.x.npc = 0.42, label.y.npc = 0.85, size = 10)
p$layers[[3]]$aes_params$textsize <- 20
p
ggsave('results/figures/fig3/brca_ssgsea_go_aging.eps', height = 9, width = 11)

