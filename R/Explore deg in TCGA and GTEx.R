source("code/useful_functions.R")
source("code/02 - survival.R")
library(pheatmap)
library(limma)
library(edgeR)
library(pbmcapply)
library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
theme_set(theme_pubr(20))

## barplot of # of DEG ##
#########################
deg <- readRDS("results/rds/logistic/deg_ensembl.rds") %>%
  dplyr::bind_rows(.id = "cancer") %>%
  dplyr::filter(cancer %in% sig_cancers) %>%
  dplyr::mutate(Sig = ifelse(adj.P.Val < 0.05, "Sig", "NS")) %>%
  dplyr::group_by(cancer) %>%
  dplyr::count(Sig) %>%
  dplyr::filter(Sig == "Sig") %>%
  dplyr::bind_rows(data.frame(cancer = setdiff(sig_cancers, .$cancer),
                              Sig = "Sig",
                              n = 0)) %>%
  dplyr::inner_join(setNames(data.frame(table(survival_data$type)), c("cancer", "samp"))) %>%
  dplyr::mutate(n = 100*n/56537,
                age_ass = ifelse(n > 1, "Age-Associated", "Not Age-Associated"))
aa_cancers <- deg$cancer[deg$age_ass == "Age-Associated"]
naa_cancers <- deg$cancer[deg$age_ass != "Age-Associated"]
total_cancers <- c(aa_cancers, naa_cancers)

ggplot(deg, aes(x = reorder(cancer, n), y = n, label = round(n,2), fill = age_ass)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Percent DEG", fill = "Molecular\nPhenotype") +
  geom_hline(yintercept = 1 , lty = 2) +
  scale_fill_manual(values = cols) +
  coord_flip() +
  theme_pubr(14) +
  theme(legend.position = c(0.8,0.2))
ggsave('results/figures/fig1/deg_counts.eps', width = 6, height = 6)

## run gsea on DEG ##
#####################
pathways.hallmark <- fgsea::gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GSEA_pathways/h.all.v7.0.symbols.gmt")
deg <- lapply(readRDS("results/rds/logistic/deg_gene.rds"), rownames_to_column, "gene")
gsea_results <- pbmclapply(deg, function(x){
  set.seed(2019)
  fgsea::fgsea(pathways = pathways.hallmark,
               stats = setNames(x %>% dplyr::arrange(t) %>% dplyr::pull(t), x %>% dplyr::arrange(t) %>% dplyr::pull(gene)),
               nperm = 1000)
}, mc.cores = 10)

gsea_results <- bind_rows(gsea_results, .id = "type") %>%
  filter(type %in% total_cancers) %>%
  mutate(cancer_anno = ifelse(type %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
         NES = ifelse(padj < 0.05, NES, 0),
         padj = ifelse(padj < 0.05, padj, NA),
         pathway = stringr::str_replace_all(substring(pathway, 10), "_", " "))
fdr <- reshape2::dcast(gsea_results, type + cancer_anno~pathway, value.var = "padj")
nes <- reshape2::dcast(gsea_results, type + cancer_anno~pathway, value.var = "NES")

rowData <- fdr[,1:2]
fdr <- fdr %>% column_to_rownames("type") %>% dplyr::select(-cancer_anno)
nes <- nes %>% column_to_rownames("type") %>% dplyr::select(-cancer_anno)

nes <- nes[, colSums(nes == 0) < 8]
nes <- nes[, setdiff(colnames(nes), c("MYC TARGETS", "G2M CHECKPOINT", "MYC TARGETS V1", "MYC TARGETS V2"))]

range <- max(abs(nes))
len <- 21
breakList <- seq(-range,range, length.out = len)

an <- HeatmapAnnotation(df = rowData %>% column_to_rownames('type') %>% dplyr::rename(Type = cancer_anno), name = 'Type', col = list(Type = cols[1:2]),
                       annotation_legend_param = list(Type = list(direction = "horizontal", ncol = 2, title_position = "topleft")))
hm <- Heatmap(t(nes), name = 'NES', rect_gp = gpar(col = "white", lwd = 2), row_names_gp = gpar(fontsize = 10), top_annotation = an, 
        heatmap_legend_param = list(title_position = "topcenter", legend_direction = "horizontal", legend_width = unit(5, "cm")),
        width = unit(14, "cm"), col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "ivory", "red")))
dev.off()
cairo_ps('results/figures/fig2/tcga_pathways_deg_not_aa_naa.eps', width = 10)
draw(hm, heatmap_legend_side = "top", annotation_legend_side = "top", merge_legend = TRUE)
dev.off()

# pheatmap(t(nes[rowData$cancer_anno == "Age-Associated", ]), scale = 'none', color = colorRampPalette(c("blue","ivory","red"))(len), breaks = breakList,
#          cellheight = 20, cellwidth = 20, fontsize = 20, treeheight_row = 15, treeheight_col = 15, 
#          border_color = "white", filename = "results/figures/fig2/tcga_pathways_deg_AA.jpeg", width = 8.5, height = 4)
# pheatmap(t(nes[rowData$cancer_anno != "Age-Associated", ]), scale = 'none', color = colorRampPalette(c("blue","ivory","red"))(len), breaks = breakList,
#          cellheight = 20, cellwidth = 20, fontsize = 20, treeheight_row = 15, treeheight_col = 15,
#          border_color = "white", filename = "results/figures/fig2/tcga_pathways_deg_notAA.jpeg", width = 10, height = 4)
keep_path <- colnames(nes); dev.off()

##  do differential expression on TCGA Normal data ##
#####################################################
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
colnames(normal_data) <- gsub('\\.', '-', substr(colnames(normal_data), 1, 12))

# get samples
brca_samples <- survival_data %>% filter(type == "BRCA", Quartiles %in% c(1,4)) %>% pull(bcr_patient_barcode)
common_samples <- intersect(colnames(normal_data), brca_samples)
surv <- survival_data %>% filter(bcr_patient_barcode %in% common_samples) %>% 
  dplyr::filter(Quartiles %in% c(1,4)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old")) %>%
  arrange(bcr_patient_barcode)
normal_data <- normal_data[, colnames(normal_data) %in% common_samples]
normal_data <- normal_data[, sort(colnames(normal_data))]
all(colnames(normal_data) == surv$bcr_patient_barcode)

# run limma
design <- model.matrix(~ 0 + surv$Quartiles)
colnames(design) <- c("Old","Young")
contmatrix <- makeContrasts("Old-Young", levels = design)
normal_data <- normal_data[which(rowSums(cpm(normal_data)> 2) >= 2) ,]
print(dim(normal_data))
Data.voom <- voom(normal_data,plot=FALSE)
fit <- lmFit(Data.voom,design)
fit2 <- contrasts.fit(fit, contmatrix)
fit2 <-eBayes(fit2)
deg.res <- topTable(fit2,n=Inf, coef = 1) %>%
  getGeneNamesFromENSELBL() %>%
  tibble::rownames_to_column("gene")

# run GSEA
set.seed(2019)
brca_normal <- fgsea::fgsea(pathways = pathways.hallmark,
                            stats = setNames(deg.res %>% dplyr::arrange(t) %>% dplyr::pull(t), deg.res %>% dplyr::arrange(t) %>% dplyr::pull(gene)),
                            nperm = 1000)


##  do differential expression on METABRIC data ##
##################################################
source("code/02 - survival.R")
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
metabric_expression <- fread("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/brca_metabric_all_Data/data_expression.txt")[,-2]
metabric_expression <- column_to_rownames(metabric_expression, "Hugo_Symbol")
colnames(metabric_expression) <- gsub('\\.', '-', colnames(metabric_expression))

# identify common cases
common_cases <- intersect(colnames(metabric_expression), metabric_clinical$X.Patient.Identifier)

# Keep common samples
metabric_expression <- metabric_expression[,common_cases]
metabric_clinical <- metabric_clinical[metabric_clinical$X.Patient.Identifier %in% common_cases, ]

# Set up for limma
###########################################################################################################
metabric_expression <- metabric_expression[, sort(colnames(metabric_expression))]

input <- data.frame(Patient = as.character(metabric_clinical$X.Patient.Identifier),
                    Age = metabric_clinical$Quartile,
                    row.names = as.character(metabric_clinical$X.Patient.Identifier))
input <- input[sort(rownames(input)), ]

# Verify this is true
all(input$Patient == colnames(metabric_expression))

# Limma
###########################################################################################################
input$Age <- factor(input$Age, levels = c("Old","Young"))
condition <- c("Old","Young")
comparison <- "Old-Young"
design <- model.matrix(~0+ input$Age)
colnames(design) <- condition
contmatrix <- makeContrasts(as.character(comparison),levels=design)
print(dim(metabric_expression))
fit <- lmFit(metabric_expression,design)
fit2 <- contrasts.fit(fit,contmatrix)
fit2 <-eBayes(fit2)
metabric.deg <- topTable(fit2,n=Inf) %>%
  rownames_to_column("gene")

set.seed(2019)
metabric.gsea <- fgsea::fgsea(pathways = pathways.hallmark,
                              stats = setNames(metabric.deg %>% dplyr::arrange(t) %>% dplyr::pull(t), 
                                               metabric.deg %>% dplyr::arrange(t) %>% dplyr::pull(gene)),
                              nperm = 1000)
dev.off()
metabric.gsea %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(., aes(x = reorder(pathway, NES), y = NES))+
  geom_bar(stat = "identity") + coord_flip() + labs(caption = "metabric")
ggsave("results/figures/fig2/metabric_gsea.eps")

# do differential expression on GTEx data
young <- c("20-29","30-39", "40-49")
old <- c("60-69", "70-79")
meta1 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t") %>%
  dplyr::filter(SMAFRZE == "RNASEQ",
                SMTS %in% c("Breast", "Thyroid", "Ovary", "Blood", "Lung", "Uterus"),
                SMTSD != "Cells - EBV-transformed lymphocytes")
meta2 <- fread("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt") %>%
  dplyr::mutate(AGE = ifelse(AGE %in% young, "Young", ifelse(AGE %in% old, "Old", NA))) %>%
  tidyr::drop_na()
gtex_mat <- data.frame(CePa::read.gct("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/full_rnaseq_matrix/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"))
colnames(gtex_mat) <- gsub("\\.", "-", colnames(gtex_mat))

# identify samples for DE
gtex_mat <- gtex_mat[, colnames(gtex_mat) %in% meta1$SAMPID]
samples_to_keep <- gsub("^([^-]*-[^-]*)-.*$", "\\1", meta1$SAMPID)
meta2 <- meta2 %>%
  dplyr::filter(SUBJID %in% samples_to_keep)
full_ids <- lapply(split(meta1 , f = meta1$SMTS), function(x) x$SAMPID)


# this is to run DE
gtex_DEG <- function(expressionMat, metaData2, complete_sample_names){
  
  expressionMat <- expressionMat[, colnames(expressionMat) %in% unlist(complete_sample_names)]
  
  # match column names with meta2
  colnames(expressionMat) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", colnames(expressionMat))
  common_cases <- intersect(metaData2$SUBJID, colnames(expressionMat))
  metaData2 <- metaData2 %>%
    dplyr::filter(SUBJID %in% common_cases) %>%
    dplyr::select(SUBJID, AGE) %>%
    dplyr::arrange(SUBJID) %>% 
    dplyr::mutate(AGE = factor(AGE, levels = c("Old", "Young")))
  expressionMat <- expressionMat[, colnames(expressionMat) %in% common_cases]
  expressionMat <- expressionMat[, sort(colnames(expressionMat))]
  
  # check
  all(colnames(expressionMat) == metaData2$SUBJID)
  
  # limma
  design <- model.matrix(~0+ metaData2$AGE)
  colnames(design) <- c("Old", "Young")
  contmatrix <- makeContrasts("Old - Young", levels = design)
  expressionMat <- expressionMat[which(rowSums(cpm(expressionMat)> 2) >= 2) ,]; print(dim(expressionMat))
  Data.voom <- voom(expressionMat, plot=FALSE)
  fit <- lmFit(Data.voom,design)
  fit2 <- contrasts.fit(fit,contmatrix)
  fit2 <-eBayes(fit2)
  tt <- topTable(fit2, n=Inf) %>%
    getGeneNamesFromENSELBL(.) %>%
    rownames_to_column("gene")
  
  return(tt)
}

# run DE
gtex_de <- pbmclapply(full_ids, gtex_DEG, expressionMat = gtex_mat, metaData2 = meta2, mc.cores = 8)

# gsea on GTEx
gtex_gsea <- pbmclapply(gtex_de, function(x){
  set.seed(2019)
  fgsea::fgsea(pathways = pathways.hallmark,
               stats = setNames(x %>% dplyr::arrange(t) %>% dplyr::pull(t), x %>% dplyr::arrange(t) %>% dplyr::pull(gene)),
               nperm = 1000)},
  mc.cores = 8)

# rerun GSEA on tcga Age Associated
deg <- lapply(readRDS("results/rds/logistic/deg_gene.rds")[aa_cancers], rownames_to_column, "gene")
gsea_results <- pbmclapply(deg, function(x){
  set.seed(2019)
  fgsea::fgsea(pathways = pathways.hallmark,
               stats = setNames(x %>% dplyr::arrange(t) %>% dplyr::pull(t), x %>% dplyr::arrange(t) %>% dplyr::pull(gene)),
               nperm = 1000)
}, mc.cores = 10)

# combine all gsea results
# use blood, laml, breast, brca, lung, lusc, ovary, ov, thca, thyroid, ucec, uterus
# sub_tissues <- c("Blood", "LAML", "Breast", "BRCA", "Lung", "LUSC", "Ovary", "OV", "Thyroid", "THCA", "UCEC", "Uterus")
sub_tissues <- c("BRCA", "Breast", "THCA", "Thyroid", "METABRIC", "BRCA.Normal")
added <- setNames(list(metabric.gsea, brca_normal), c("METABRIC", "BRCA.Normal"))
combined_gsea <- rbind(bind_rows(bind_rows(gtex_gsea, .id = "tissue"), bind_rows(gsea_results, .id = "tissue"), bind_rows(added, .id = "tissue"))) %>%
  dplyr::mutate(Status = ifelse(tissue %in% c("BRCA", "THCA", "METABRIC"), "Tumor", "Normal"),
                NES = ifelse(padj < 0.05, NES, 0),
                padj = ifelse(padj < 0.05, padj, NA),
                pathway = stringr::str_replace_all(substring(pathway, 10), "_", " ")) %>%
  dplyr::filter(tissue %in% sub_tissues)
fdr <- reshape2::dcast(combined_gsea, tissue + Status~pathway, value.var = "padj")
nes <- reshape2::dcast(combined_gsea, tissue + Status~pathway, value.var = "NES")

rowData <- fdr[,1:2] %>%
  column_to_rownames("tissue")
fdr <- fdr %>% column_to_rownames("tissue") %>% dplyr::select(-Status)
nes <- nes %>% column_to_rownames("tissue") %>% dplyr::select(-Status)

# nes <- nes[, colSums(nes == 0) < 10]
nes <- nes[, colnames(nes) %in% keep_path]

an = HeatmapAnnotation(df = anno, name = 'Source', col = list(Source = cols[5:6]),
                       annotation_legend_param = list(Source = list(direction = "horizontal", ncol = 2, title_position = "topleft")))
hm = Heatmap(t(nes), name = 'NES', top_annotation = an, 
             col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "ivory", "red")), rect_gp = gpar(col = "white", lwd = 2), row_names_gp = gpar(fontsize = 10),
             heatmap_legend_param = list(title_position = "topcenter", legend_width = unit(5, "cm"), legend_direction = 'horizontal',
                                         at = c(-4, -2, 0, 2, 4), labels = c('Young', -2, 0, 2, 'Old')), 
             width = unit(7, 'cm'))

cairo_ps('results/figures/fig2/cancer_vs_normal_gsea.eps', width = 6)
draw(hm, heatmap_legend_side = 'top', merge_legend = T)
dev.off()

# anno = data.frame(Source = c('Tumor', 'Normal', 'Normal', 'Tumor', 'Tumor', 'Normal'), row.names = colnames(t(nes)))
# 
# pheatmap(t(nes), scale = 'none', annotation_col = rowData, color = colorRampPalette(c("blue","ivory","red"))(len), breaks = breakList,
#          cellheight = 20, cellwidth = 20, fontsize = 20, treeheight_row = 15, treeheight_col = 15,
#          border_color = "white", annotation_colors = list(Status = c(`Tumor`="black", `Normal` = "gray")),
#          width = 10, height = 5, cutree_cols = 2)
# 
# dev.off()
# pheatmap(t(nes), scale = 'none', annotation_col = rowData, color = colorRampPalette(c("blue","ivory","red"))(len), breaks = breakList,
#          cellheight = 20, cellwidth = 20, fontsize = 20, treeheight_row = 15, treeheight_col = 15,
#          border_color = "white", annotation_colors = list(Status = c(`Tumor`="black", `Normal` = "gray")),
#          width = 10, height = 5, cutree_cols = 2, filename = "results/figures/fig2/cancer_vs_normal_gsea.jpeg")
