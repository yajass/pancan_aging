source("code/useful_functions.R")
source("code/02 - survival.R")
library(TCGAutils); library(TCGAmutations); library(TCGAbiolinks)
library(maftools)
library(dplyr); library(stringr)
library(ggsci); library(ggpubr)

# Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
          "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
          "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")

# load and merge mafs
dir <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/tcga_mutations/"
aa_maf <- list.files(dir); aa_maf <- aa_maf[gsub('.{6}$', '', aa_maf) %in% aa_cancers]
lapply(paste0(dir, aa_maf), load, temp_env <- new.env())
env_list <- as.list.environment(temp_env)
aa_maf <- merge_mafs(env_list)
naa_maf <- list.files(dir); naa_maf <- naa_maf[gsub('.{6}$', '', naa_maf) %in% naa_cancers]
lapply(paste0(dir, naa_maf), load, temp_env <- new.env())
env_list <- as.list.environment(temp_env)
naa_maf <- merge_mafs(env_list)
rm(temp_env); rm(env_list)

# load METABRIC maf
young_cutoff = survival_data %>% filter(type == "BRCA", Quartiles == 1) %>% pull(age_at_initial_pathologic_diagnosis) %>% max()
old_cutoff = survival_data %>% filter(type == "BRCA", Quartiles == 4) %>% pull(age_at_initial_pathologic_diagnosis) %>% min()
clinical <- fread("../Data/brca_metabric_all_Data/data_clinical_patient copy.txt", skip = 2) %>%
  dplyr::mutate(Quartiles = ifelse(AGE_AT_DIAGNOSIS <= young_cutoff, "Young",
                                   ifelse(AGE_AT_DIAGNOSIS >= old_cutoff, "Old", "Middle.Aged"))) %>%
  dplyr::rename(Tumor_Sample_Barcode = PATIENT_ID) %>%
  dplyr::select(Tumor_Sample_Barcode, Quartiles) %>%
  mutate(type = "BRCA", Dataset = "METABRIC")
metabric_maf = read.maf("../Data/brca_metabric_all_Data/data_mutations_extended.txt",
                        isTCGA = F, clinicalData = clinical)
metabric_maf = subsetMaf(metabric_maf, clinQuery = "Quartiles %in% c('Young', 'Old')")

# merge with quartile information
aa_maf@clinical.data <- aa_maf@clinical.data %>% 
  inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles, type) %>% dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", NA)),
         Dataset = "TCGA") %>%
  filter(Quartiles %in% c("Young", "Old")) %>%
  mutate(Quartiles = factor(Quartiles, levels = c("Young", "Old"))) %>%
  dplyr::select(Tumor_Sample_Barcode, Dataset, type, Quartiles)
aa_maf <- subsetMaf(aa_maf, clinQuery = "Quartiles %in% c('Young', 'Old')")
try.aa_maf <- merge_mafs(list(TCGA=aa_maf, METABRIC=metabric_maf), verbose = T)
naa_maf@clinical.data <- naa_maf@clinical.data %>% 
  inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles, type) %>% dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", NA)),
         Dataset = "TCGA") %>%
  filter(Quartiles %in% c("Young", "Old")) %>%
  mutate(Quartiles = factor(Quartiles, levels = c("Young", "Old"))) %>%
  dplyr::select(Tumor_Sample_Barcode, Dataset, type, Quartiles)
naa_maf <- subsetMaf(naa_maf, clinQuery = "Quartiles %in% c('Young', 'Old')")
saveRDS(aa_maf, "results/rds/aa_maf.rds")
saveRDS(naa_maf, "results/rds/naa_maf.rds")
saveRDS(try.aa_maf, 'results/rds/try.aa_maf.rds')

# set colors
age_col <- setNames(c(pal_jco()(2)[1], pal_jco()(2)[2]), c("Young", "Old"))
cancer_col_aa <- setNames(RColorBrewer::brewer.pal(n = length(unique(aa_maf@clinical.data$type)), name = 'Dark2'), unique(aa_maf@clinical.data$type))
cancer_col_naa <- setNames(pal_npg()(length(unique(naa_maf@clinical.data$type))), unique(naa_maf@clinical.data$type))
col_dataset <- setNames(pal_jama()(2), c("TCGA", "METABRIC"))
colnames(aa_maf@clinical.data)[4] <- "Age"; colnames(naa_maf@clinical.data)[4] <- "Age"
aa_maf@clinical.data$Age <- factor(aa_maf@clinical.data$Age, levels = c("Young", "Old"))
naa_maf@clinical.data$Age <- factor(naa_maf@clinical.data$Age, levels = c("Young", "Old"))

tiff("results/figures/fig4/oncoplot_aa.tiff", res = 320, units = "in", width = 8, height = 7)
oncoplot(maf = aa_maf, 
         top = 17, genesToIgnore = flags, 
         # bgCol = "ivory",
         legendFontSize = 1,sepwd_samples = 0,draw_titv = F,logColBar = T,
         sortByAnnotation = T,annotationFontSize = 1, fontSize = 1, titleFontSize = 1, 
         annotationColor = list(Age = age_col, type = cancer_col_aa),drawRowBar = T,
         gene_mar = 7, showTitle = FALSE, # writeMatrix = TRUE,
         clinicalFeatures = c("Age","type")) # which columns should be plotted
dev.off()

tiff("results/figures/fig4/oncoplot_naa.tiff", res = 320, units = "in", width = 8, height = 7)
oncoplot(maf = naa_maf, 
         top = 14, genesToIgnore = flags, 
         # bgCol = "white",
         legendFontSize = 1,sepwd_samples = 0,draw_titv = F,logColBar = T,
         sortByAnnotation = T,annotationFontSize = 1, fontSize = 1, titleFontSize = 1,
         annotationColor = list(Age = age_col, type = cancer_col_naa),drawRowBar = T,
         gene_mar = 7, showTitle = FALSE, # writeMatrix = TRUE,
         clinicalFeatures = c("Age","type")) # which columns should be plotted
dev.off()




multiHit <- function(matrixDir, maf){
  mat = fread(matrixDir)
  mat = melt(mat, "Gene") %>%
    dplyr::rename(Tumor_Sample_Barcode = variable) %>%
    inner_join(maf@clinical.data)
  mat = mat %>%
    group_by(Gene) %>%
    dplyr::count(Age, value=="Multi_Hit") %>%
    drop_na %>% dplyr::rename(multihit = `value == "Multi_Hit"`) %>% ungroup()
  mat = split(mat, mat$Gene)
  mat = mat[(lapply(mat, nrow) %>% unlist) != 2]
  
  res = lapply(mat, function(x) {
    m = x %>% dplyr::select(-Gene) %>% 
      reshape2::dcast(multihit~Age, value.var="n") %>%
      column_to_rownames("multihit")
    m[is.na(m)] = 0
    return(rstatix::fisher_test(m, detailed = TRUE))
  }) %>% bind_rows() %>% mutate(Gene = names(mat)) %>%
    rstatix::adjust_pvalue("p","fdr") %>%
    rstatix::add_significance(p.col = "fdr") %>%
    filter(fdr < 0.05)
  
  plt = mat[res$Gene] %>% bind_rows() %>%
    mutate(multihit = ifelse(multihit, "Multi-Hit", "Single Hit")) %>%
    group_by(Gene,Age) %>%
    mutate(n = 100*n/sum(n)) %>%
    filter(multihit == "Multi-Hit")
  return(plt)
}

multiHit("results/spreadsheets/oncoprint/aa_onco_matrix.txt", aa_maf) %>%
  ggplot(., aes(x = Gene, y = n, fill=Age)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  theme_pubr() + labs(x = NULL, y = "% Multi-Hit") +
  scale_fill_jco() + coord_flip() +
  theme(axis.ticks.x = element_blank(),
        legend.position = 'none')
ggsave("results/figures/fig4/aa_multihit.eps", height = 6, width = 2)


multiHit("results/spreadsheets/oncoprint/naa_onco_matrix.txt", naa_maf) %>%
  ggplot(., aes(x = Gene, y = n, fill=Age)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  theme_pubr(20) + labs(x = NULL, y = "% Multi-Hit") +
  scale_fill_jco() + coord_flip() +
  theme(axis.ticks.x = element_blank(),
        legend.position = 'none')
# ggsave("results/figures/fig4/naa_multihit.eps")  ## no sig results




# source("code/useful_functions.R")
# source("code/02 - survival.R")
# library(TCGAutils); library(TCGAmutations); library(TCGAbiolinks)
# library(maftools)
# library(dplyr); library(stringr)
# library(ggsci); library(ggpubr)
# 
# # Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
# flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
#           "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
#           "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")
# 
# # load and merge mafs
# dir <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/tcga_mutations/"
# aa_maf <- list.files(dir); aa_maf <- aa_maf[gsub('.{6}$', '', aa_maf) %in% aa_cancers]
# lapply(paste0(dir, aa_maf), load, temp_env <- new.env())
# env_list <- as.list.environment(temp_env)
# aa_maf <- merge_mafs(env_list)
# naa_maf <- list.files(dir); naa_maf <- naa_maf[gsub('.{6}$', '', naa_maf) %in% naa_cancers]
# lapply(paste0(dir, naa_maf), load, temp_env <- new.env())
# env_list <- as.list.environment(temp_env)
# naa_maf <- merge_mafs(env_list)
# rm(temp_env); rm(env_list)
# 
# # merge with quartile information
# aa_maf@clinical.data <- aa_maf@clinical.data %>% 
#   inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles, type) %>% dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode)) %>%
#   mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", NA))) %>%
#   filter(Quartiles %in% c("Young", "Old")) %>%
#   mutate(Quartiles = factor(Quartiles, levels = c("Young", "Old")))
# aa_maf <- subsetMaf(aa_maf, clinQuery = "Quartiles %in% c('Young', 'Old')")
# naa_maf@clinical.data <- naa_maf@clinical.data %>% 
#   inner_join(survival_data %>% dplyr::select(bcr_patient_barcode, Quartiles, type) %>% dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode)) %>%
#   mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", NA))) %>%
#   filter(Quartiles %in% c("Young", "Old")) %>%
#   mutate(Quartiles = factor(Quartiles, levels = c("Young", "Old")))
# naa_maf <- subsetMaf(naa_maf, clinQuery = "Quartiles %in% c('Young', 'Old')")
# saveRDS(aa_maf, "results/rds/aa_maf.rds")
# saveRDS(naa_maf, "results/rds/naa_maf.rds")
# 
# # set colors
# age_col <- setNames(c(pal_jco()(2)[1], pal_jco()(2)[2]), c("Young", "Old"))
# cancer_col_aa <- setNames(RColorBrewer::brewer.pal(n = length(unique(aa_maf@clinical.data$type)), name = 'Dark2'), unique(aa_maf@clinical.data$type))
# cancer_col_naa <- setNames(pal_npg()(length(unique(naa_maf@clinical.data$type))), unique(naa_maf@clinical.data$type))
# colnames(aa_maf@clinical.data)[228] <- "Age"; colnames(naa_maf@clinical.data)[228] <- "Age"
# aa_maf@clinical.data$Age <- factor(aa_maf@clinical.data$Age, levels = c("Old", "Young"))
# 
# oncoplot(maf = aa_maf, 
#          top = 15, genesToIgnore = flags, 
#          # bgCol = "white",
#          legendFontSize = 1,sepwd_samples = 0,draw_titv = F,logColBar = T,
#          sortByAnnotation = T,annotationFontSize = 1, fontSize = 1, titleFontSize = 1,
#          annotationColor = list(Age = age_col, type = cancer_col_aa),
#          gene_mar = 7, showTitle = FALSE,
#          clinicalFeatures = c("Age","type")) # which columns should be plotted
# 
# # multi-hit code from github
# xx <- aa_maf@data %>%
#   dplyr::filter(Hugo_Symbol %in% c('TP53','PIK3CA','RYR2','PTEN','CSMD3','IDH1','ZFHX4','LRP1B','ARID1A','ATRX')) %>%
#   dplyr::group_by(Tumor_Sample_Barcode) %>%
#   inner_join(aa_maf@clinical.data %>% dplyr::select(Tumor_Sample_Barcode, Age)) %>%
#   dplyr::count(Tumor_Sample_Barcode, Age, Hugo_Symbol)
# 
# pvals <- numeric(length = length(unique(xx$Hugo_Symbol)))
# f.t <- list()
# res <- list()
# for (i in 1:length(unique(xx$Hugo_Symbol))) {
#   df <- xx %>% dplyr::filter(Hugo_Symbol == unique(xx$Hugo_Symbol)[i])
#   tab <- setNames(data.frame(matrix(c(sum(df$n[df$Age == "Young"] == 1),
#                                       sum(df$n[df$Age == "Old"] == 1),
#                                       sum(df$n[df$Age == "Young"] > 1),
#                                       sum(df$n[df$Age == "Old"] > 1)),
#                                     nrow = 2, ncol = 2),
#                              row.names = c("Young", "Old")),
#                   c("Single", "Multi"))
#   print(unique(xx$Hugo_Symbol)[i])
#   res[[i]] <- tab
#   f.t[[i]] <- fisher.test(tab)
#   pvals[i]<- fisher.test(t(tab))$p.value
#   names(pvals)[i] <- unique(xx$Hugo_Symbol)[i]
#   names(f.t)[[i]] <- unique(xx$Hugo_Symbol)[i]
#   names(res)[[i]] <- unique(xx$Hugo_Symbol)[i]
# }
# 
# res <- lapply(res, function(x){
#   x %>% rownames_to_column("age") %>% melt("age")
# }) %>% bind_rows() %>% mutate(gene = rep(names(res), each=4))
# 
# res %>%
#   filter(gene %in% names(p.adjust(pvals, "fdr") < 0.05)[p.adjust(pvals, "fdr") < 0.05]) %>%
#   group_by(gene, variable) %>%
#   mutate(perc = 100*value/sum(value),
#          age = factor(age, levels = c("Young", "Old")),
#          variable = ifelse(variable == "Single", "Single\nHit", "Multi\nHit")) %>%
#   ggbarplot(., x = "variable", y = "perc", fill = "age", label = TRUE, facet.by = "gene", lab.nb.digits = 2, lab.pos = "in", palette = "jco",lab.size = 5) +
#   labs(x = NULL, y = "Percent of Whole", fill = NULL) + theme_pubr(20, legend = "none")
