options(stringsAsFactors = FALSE)

# Load everything
################################################################
# load rna age prection model
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/Forked_Packages/RNAAgeCalc-master/R/predict_age.R")
source("/Users/Yajas/Documents/Elemento/tcga_aging_final/Forked_Packages/RNAAgeCalc-master/R/makeplot.R")
load("/Users/Yajas/Documents/Elemento/tcga_aging_final/Forked_Packages/RNAAgeCalc-master/R/sysdata.rda")
source("code/useful_functions.R")
source("code/02 - survival.R")

# load packages
library(impute)
library(pbmcapply)
library(rstatix)
library(patchwork)
library(ggpubr); theme_set(theme_pubr())
library(tidyverse)
library(stringr)

# Create sample table
################################################################
samples <- data.frame(name = c(aa_cancers, naa_cancers),
                      tissue = c("breast","brain","lung","ovary","thyroid","uterus",
                                 "bladder","colon","brain","salivary","kidney","blood","pancreas","kidney","adipose_tissue", "stomach"))
fpkm_dir <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/schults_fpkm_combat/"
files <- list.files(fpkm_dir)
t.files <- files[endsWith(files, "-t.txt")]
n.files <- files[endsWith(files, "tcga.txt")]

# match schultz data and expand table
samples$gtex_fpkm <- lapply(samples$tissue, function(x) paste0(fpkm_dir, files[str_detect(string = files, pattern = x)])) %>% unlist()
samples$tcga_t_fpkm <- lapply(tolower(samples$name), function(x) paste0(fpkm_dir, t.files[str_detect(string = t.files, pattern = x)])) %>% unlist()
samples$tcga_n_fpkm <- lapply(tolower(samples$name), function(x) paste0(fpkm_dir, n.files[str_detect(string = n.files, pattern = x)])) %>% unlist()
samples <- samples[apply(samples[,3:5], 2, endsWith, ".txt")[,1], ]
samples$tissue[samples$tissue == "salivary"] <- "salivary_gland"
# samples <- samples[-9, ]
samples <- data.table(samples); setkey(samples, 'name')


# Functions for age prediction
################################################################
calculateAgeTCGA <- function(tumor_type, metadata = survival_data, sig = "Dev", stype = "all"){
  
  # load fpkm
  tcga_t <- fread(samples[tumor_type]$tcga_t_fpkm)[,-2] %>% column_to_rownames("Hugo_Symbol")
  colnames(tcga_t) <- gsub('\\.','-',substring(colnames(tcga_t),1,12))
  common_t <- intersect(colnames(tcga_t), survival_data$bcr_patient_barcode)
  tcga_t <- tcga_t[,common_t]
  tcga_t <- tcga_t[, sort(colnames(tcga_t))]
  tcga_n <- fread(samples[tumor_type]$tcga_n_fpkm)[,-2] %>% column_to_rownames("Hugo_Symbol")
  colnames(tcga_n) <- gsub('\\.','-',substring(colnames(tcga_n),1,12))
  common_n <- intersect(colnames(tcga_n), survival_data$bcr_patient_barcode)
  tcga_n <- tcga_n[,common_n]
  tcga_n <- tcga_n[, sort(colnames(tcga_n))]
  
  # get metadata
  # match1 <- match(survival_data$bcr_patient_barcode, colnames(tcga_t))
  metadata_t <- survival_data %>% dplyr::filter(bcr_patient_barcode %in% common_t) %>%
    dplyr::select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis) %>%
    dplyr::rename(UID = bcr_patient_barcode, chonage = age_at_initial_pathologic_diagnosis) %>% dplyr::arrange(UID)
  metadata_n <- survival_data %>% dplyr::filter(bcr_patient_barcode %in% common_n) %>%
    dplyr::select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis) %>%
    dplyr::rename(UID = bcr_patient_barcode, chonage = age_at_initial_pathologic_diagnosis) %>% dplyr::arrange(UID)
  
  # make sure everything matches up
  all(colnames(tcga_t) == metadata_t$UID); all(colnames(tcga_n) == metadata_n$UID)
  
  # predict age
  pred_t <- predict_age(exprdata = tcga_t, tissue = samples[tumor_type]$tissue, chronage = metadata_t,
                        exprtype = "FPKM", idtype = "SYMBOL", signature = sig, stype = stype) %>%
    rownames_to_column("SampleID") %>%
    dplyr::mutate(type = tumor_type,
                  dataset = "TCGA",
                  sample = "Tumor")
  pred_n <- predict_age(exprdata = tcga_n, tissue = samples[tumor_type]$tissue, chronage = metadata_n,
                        exprtype = "FPKM", idtype = "SYMBOL", signature = sig, stype = stype) %>%
    rownames_to_column("SampleID") %>%
    dplyr::mutate(type = tumor_type,
                  dataset = "TCGA",
                  sample = "Normal")
  
  return(bind_rows(pred_t, pred_n))
  
}

calculateAgeGTEx <- function(tumor_type, sig = "Dev", stype=  "all"){
  
  # load fpkm
  gtex <- fread(samples[tumor_type]$gtex_fpkm)[,-2] %>% column_to_rownames("Hugo_Symbol")
  colnames(gtex) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", colnames(gtex))
  
  # predict
  pred_gtex <- predict_age(exprdata = gtex, tissue = samples[tumor_type]$tissue,
                           exprtype = "FPKM", idtype = "SYMBOL", signature = sig, stype) %>%
    rownames_to_column("SampleID") %>%
    dplyr::mutate(type = tumor_type,
                  dataset = "GTEx",
                  sample = "Normal")
  return(pred_gtex)
  
}


# Predict age
################################################################
sig = "Dev" # options are "DESeq"       "Pearson"     "Dev"         "deMagalhaes" "GenAge"      "GTExAge"     "Peters"      "all"
tcga_age <- pbmclapply(samples$name, calculateAgeTCGA, metadata = survival_data, stype = 'all', mc.cores = 10) %>% bind_rows()
gtex_age <- pbmclapply(samples$name, calculateAgeGTEx, stype = 'all', mc.cores = 10) %>% bind_rows()

# Analyze results
################################################################
paired_samps <- tcga_age %>% dplyr::count(SampleID) %>% filter(n>1) %>% pull(SampleID)
# tcga_age %>% filter(SampleID %in% paired_samps) %>%
#   mutate(AgeAccDiff = RNAAge - ChronAge) %>%
#   dplyr::select(SampleID, RNAAge, ChronAge, AgeAccDiff, AgeAccelResid, type, dataset, sample) %>%
#   fwrite('results/tables/rna_age.tsv', sep = '\t', row.names = F, quote = F)

# Scatterplots for GTEx and TCGA
meta2 <- read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                  stringsAsFactors = F, header = T, sep = "\t") %>% 
  mutate(SUBJID = gsub('\\-','.',SUBJID), AGE = factor(AGE)) %>% dplyr::rename(SampleID = SUBJID)
df = data.frame(type = c('BLCA', 'BRCA', 'COAD', 'HNSC', 'KICH', 'KIRC', 'LUSC', 'THCA', 'UCEC'),
                organ = c('Bladder', 'Breast', 'Colon', 'Salivary Gland', 'Kidney', 'Kidney', 'Lung', 'Thyroid', 'Uterus'),
                stringsAsFactors = F) %>%
  mutate(organ = factor(organ, levels = c('Bladder', 'Breast', 'Colon', 'Salivary Gland', 'Kidney', 'Lung', 'Stomach', 'Thyroid', 'Uterus')))
gtex = inner_join(gtex_age, meta2) %>% inner_join(df)

tcga.plt <- ggscatter(tcga_age, x = 'ChronAge', y = 'RNAAge', alpha = 2/3) + facet_grid(vars(sample), vars(type)) + stat_cor(label.sep = '\n') + 
  labs(x = 'Age') + theme_bw(20) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
gtex.plt <- ggscatter(gtex, x = 'AGE', y = 'RNAAge', alpha = 2/3) + facet_grid(vars(dataset), vars(organ)) + labs(x = "Age") + stat_cor(method = 'spearman') +
  theme_bw(20) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

tcga.plt; ggsave('results/figures/fig3/validation_rna_age_tcga.jpeg', width = 15, height = 8, dpi = 300)
gtex.plt; ggsave('results/figures/fig3/validation_rna_age_gtex.jpeg', width = 15, height = 4.3, dpi = 300)


## Write tables
tab = rbind(gtex %>% dplyr::select(SampleID, dataset, type, sample, organ, AGE, RNAAge) %>% mutate(AGE = as.character(AGE)),
      tcga_age %>% mutate(organ = NA, ChronAge = as.character(ChronAge)) %>% dplyr::select(SampleID, dataset, type, sample, organ, ChronAge, RNAAge) %>% dplyr::rename(AGE = ChronAge))
fwrite(tab, 'results/tables/rna_age.tsv', sep = '\t', row.names = F, quote = F)

# Scaled age acceleration for paired samples
cols = setNames(ggsci::pal_aaas()(2), c('Tumor', 'Normal'))
tcga_age %>%
  dplyr::filter(SampleID %in% paired_samps) %>%
  mutate(AgeAccDiff = RNAAge - ChronAge,
         class = factor(ifelse(type %in% aa_cancers, "AA", "NAA"), levels = c('NAA', 'AA'))) %>%
  group_by(type) %>%
  mutate(AgeAccDiff = scale(AgeAccDiff)) %>%
  arrange(sample) %>%
  ggscatter(., x = "ChronAge", y = "AgeAccDiff", color = "sample", 
            add = "reg.line", conf.int = T, alpha = 1/3) +
  stat_cor(aes(color = sample), show.legend = FALSE) + scale_color_manual(values = cols) +
  labs(x = "Age", y = "Scaled RNA Age Acc Diff", color = NULL) +
  theme_pubr(20, legend = 'none') + guides(fill = F) +
  facet_wrap(~class, scales = "free_x")
ggsave("results/figures/fig3/scaled_rna_age_tn.jpeg", width = 8, height = 5, dpi = 320)


# Compare AA vs NAA age acceleration
st = tcga_age %>%
  inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles) %>% 
               dplyr::rename(SampleID = bcr_patient_barcode) %>%
               mutate(SampleID = gsub('-','.',SampleID))) %>%
  mutate(AgeAccDiff = RNAAge - ChronAge,
         class = ifelse(type %in% aa_cancers, "AA", "NAA")) %>%
  group_by(type) %>%
  mutate(AgeAccDiff = scale(AgeAccDiff)) %>%
  filter(Quartiles %in% c(1, 4), sample == 'Tumor') %>%
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', 'Old'), c('Young','Old'))) %>% ungroup() %>% group_by(Quartiles) %>%
  rstatix::wilcox_test(AgeAccDiff ~ class) %>% add_significance()
plt1 = tcga_age %>%
  inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles) %>% 
               dplyr::rename(SampleID = bcr_patient_barcode) %>%
               mutate(SampleID = gsub('-','.',SampleID))) %>%
  mutate(AgeAccDiff = RNAAge - ChronAge,
         class = ifelse(type %in% aa_cancers, "AA", "NAA")) %>%
  group_by(type) %>%
  mutate(AgeAccDiff = scale(AgeAccDiff)) %>%
  filter(Quartiles %in% c(1, 4), sample == 'Tumor') %>%
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', 'Old'), c('Young','Old'))) %>%
  ggviolin(x = 'class', y = 'AgeAccDiff', facet.by = 'Quartiles', add = 'boxplot', fill = 'class') + 
  scale_fill_brewer(palette = 'Dark2') + theme_pubr(25, legend = 'none') +
  labs(x = NULL, y = 'Scaled RNA AgeAccDiff') +
  stat_pvalue_manual(data = st, y.position = 2.5, size = 10, bracket.size =  0)

# Compare tumor vs normal age acceleration faceted by QUartile and AA status
plt2 = tcga_age %>%
  dplyr::filter(SampleID %in% paired_samps) %>%
  mutate(class = factor(ifelse(type %in% aa_cancers, "AA", "NAA"), levels = c('NAA', 'AA')),
         SampleID = gsub('\\.','-', SampleID)) %>%
  inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles) %>% dplyr::rename(SampleID = bcr_patient_barcode)) %>%
  group_by(type) %>%
  mutate(AgeAcc = scale(RNAAge - ChronAge),
         sample = factor(sample, levels = c('Normal', 'Tumor'))) %>%
  filter(Quartiles %in% c(1, 4)) %>%
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', 'Old'), levels = c('Young', 'Old'))) %>%
  ggpaired(x = 'sample', y = 'AgeAcc', id = 'SampleID', facet.by = c('class', 'Quartiles'), fill = 'sample', line.color = scales::alpha("gray", 2/3)) + 
  stat_compare_means(method = 'wilcox', label = 'p.signif', hide.ns = T, label.y = 2, label.x = 1.43, size = 10) +
  labs(x = NULL, y = 'Scaled RNA Age Acc Diff') + scale_fill_manual(values = cols) + theme_pubr(25, legend = 'none')

plt1 + plt2 + plot_layout(widths = c(0.8, 1))
ggsave('results/figures/fig3/rna_age_arranged.jpeg', width = 14, height = 5.5, dpi = 320)



tcga_age %>%
  ggscatter(x = 'ChronAge', y = 'RNAAge', facet.by = 'type', color = 'sample', scales = 'free', add = 'reg.line', palette = 'aaas', alpha = 1/2) +
  stat_cor(aes(color = sample), show.legend = F) + theme_pubr(10) + labs(x = 'Age', y = 'RNA Age', color = NULL)
ggsave("results/figures/fig3/rna_age_scatter_by_type.jpg", width = 12, height = 12)

