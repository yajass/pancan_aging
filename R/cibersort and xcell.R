library(ggsci)
library(pbmcapply)
library(data.table)
library(tidymodels)
library(TCGAutils)
source("code/useful_functions.R")
source("code/02 - survival.R")
# theme_set(theme_bw(20))

total_cancers <- c(aa_cancers, naa_cancers)
survival_data <- survival_data %>%
  dplyr::filter(Quartiles %in% c(1,4)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))

## CIBERSORT ABSOLUTE ##
########################
# read data
readData <- function(cancer, clinical, enrichment){
  # read data
  folder <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA_CIBERSORT_ABS_WCM/"
  abs <- fread(paste0(folder, "TCGA-", cancer, "_FPKM_CIBERSORT_ABS.txt"), check.names = TRUE) %>%
    dplyr::select(-P.value, -Correlation, -RMSE, -Absolute.score..sig.score.) # can try filtering by p value
  abs <- abs[TCGAbiospec(abs$TCGA_barcode)$sample %in% c('01','03'), ] %>%
    dplyr::mutate(TCGA_barcode = TCGAbiospec(TCGA_barcode)$submitter_id)
  abs$Lymphocyte.score <- rowSums(abs[,2:13])
  abs$Leukocyte.score <- rowSums(abs[,2:ncol(abs)])
  cell_types <- colnames(abs)[-1]
  full_df <- merge(abs, clinical, by.x = "TCGA_barcode", by.y = "bcr_patient_barcode")
  return(list(df = full_df, cell_types=cell_types))
}

# run a linear model for cell ~ age + confounders, run pearson correlation
cibersortlm <- function(cell_type, complete_df, conf){
  # setup confounder formula
  conf = conf[conf$type == unique(complete_df$type),]$logistic
  if(conf == ""){
    conf = NULL
    form1 <- as.formula("response ~  Quartiles")
  } else {
    form1 <- as.formula(paste0("response ~ Quartiles +", conf))
  }
  
  # get resposne
  response <- complete_df %>% dplyr::select(cell_type) %>% pull()
  
  # make response normally distributed
  response <- RNOmni::rankNorm(response)
  
  complete_df$Quartiles <- factor(complete_df$Quartiles, levels = c("Young", "Old"))
  # make lm
  mod <- lm(form1, data = complete_df)
  
  
  # do correlation
  corr <- cor.test(response, ifelse(complete_df$Quartiles == "Young", 0, 1), method = "spearman")
  cor.df <- data.frame(method = "Spearman", cor = corr$estimate, p = corr$p.value, 
                       cancer = unique(complete_df$type), cell_type = cell_type,
                       row.names = NULL)
  return(list(lm_res = tidy(mod) %>% 
                dplyr::mutate(cancer = unique(complete_df$type),
                              cell_type = cell_type),
              cor_res = cor.df))
}

# run analysis for total cancers
abs.data <- setNames(lapply(total_cancers, readData, clinical = survival_data), total_cancers)
for (i in seq_along(total_cancers)) {
  res <- pbmclapply(abs.data[[i]]$cell_types, cibersortlm, abs.data[[i]]$df, conf = confounders, mc.cores = 10)
  if(i == 1){
    lm.res = lapply(res, pluck, 1) %>% bind_rows()
    cor.res = lapply(res, pluck, 2) %>% bind_rows()
  } else{
    lm.res = bind_rows(lm.res, lapply(res, pluck, 1) %>% bind_rows())
    cor.res = bind_rows(cor.res, lapply(res, pluck, 2) %>% bind_rows())
  }
}

cibersort.abs.lm <- lm.res %>%
  filter(term == "QuartilesOld",
         cell_type == "Lymphocyte.score") %>%
  rstatix::adjust_pvalue(p.col = "p.value", method = "fdr")
cibersort.abs.corr <- cor.res %>%
  filter(cell_type == "Lymphocyte.score") %>%
  rstatix::adjust_pvalue(p.col = "p", method = "fdr")

## run xCell ##
###############
## run for all cancers
# get data
xcell <- fread("https://xcell.ucsf.edu/xCell_TCGA_RSEM.txt", check.names = TRUE) %>%
  column_to_rownames("V1") %>% t(.) %>% as_tibble(xcell, rownames = "TCGA_barcode") %>%
  dplyr::mutate(TCGA_barcode = as.character(gsub("\\.","-", TCGA_barcode)))
xcell <- xcell[endsWith(xcell$TCGA_barcode, '01') |endsWith(xcell$TCGA_barcode, '03'), ] %>%
  dplyr::mutate(TCGA_barcode = (substring(TCGA_barcode, 1,12)))
cell_types <- colnames(xcell)[-1]
combined_df <- merge(xcell, survival_data, by.x = "TCGA_barcode", by.y = "bcr_patient_barcode", all.x = TRUE) %>%
  dplyr::filter(type %in% total_cancers) %>%
  dplyr::mutate(type = factor(type, levels = total_cancers))
combined_df <- split(combined_df, f = combined_df$type)
all(as.character(unname(unlist(lapply(combined_df, function(x) unique(x$type))))) == total_cancers)

# run analysis
for (i in seq_along(total_cancers)) {
  res <- pbmclapply(cell_types, cibersortlm, combined_df[[i]], conf = confounders, mc.cores = 10)
  if(i == 1){
    lm.res = lapply(res, pluck, 1) %>% bind_rows()
    cor.res = lapply(res, pluck, 2) %>% bind_rows()
  } else{
    lm.res = bind_rows(lm.res, lapply(res, pluck, 1) %>% bind_rows())
    cor.res = bind_rows(cor.res, lapply(res, pluck, 2) %>% bind_rows())
  }
}

xcell.lm <- 
  lm.res %>%
  dplyr::filter(term == "QuartilesOld") %>%
  dplyr::group_by(cancer) %>%
  dplyr::mutate(fdr = p.adjust(p.value, "fdr"),
                dir = ifelse(estimate < 0 & fdr < 0.05, "Younger",
                             ifelse(estimate > 0 & fdr < 0.05, "Older",
                                    "NS")),
                dir = factor(dir, levels = c("Younger", "Older", "NS")),
                variable = paste0(cancer,"-",cell_type),
                cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  dplyr::filter(fdr < 0.05,
                cell_type %in% cell_types[c(1,4,6:15,17,20,22,29,32:34,38,41,46,48,49,50,57,61:64)])

xcell.corr <- cor.res %>%
  dplyr::group_by(cancer) %>%
  dplyr::mutate(fdr = p.adjust(p, "fdr"),
                dir = ifelse(cor < 0 & fdr < 0.05, "Younger",
                             ifelse(cor > 0 & fdr < 0.05, "Older",
                                    "NS")),
                dir = factor(dir, levels = c("Younger", "Older", "NS")),
                variable = paste0(cancer,"-",cell_type),
                cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  dplyr::filter(fdr < 0.05,
                cell_type %in% cell_types[c(1,4,6:15,17,20,22,29,32:34,38,41,46,48,49,50,57,61:64)])

## CIBERSORT RELATIVE ##
########################
# load data
rel <- fread("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Immune Landscape of Cancer/TCGA.Kallisto.fullIDs.cibersort.relative.tsv")[,-2] %>%
  dplyr::rename(TCGA_barcode = SampleID) %>%
  dplyr::select(-P.value, -Correlation, -RMSE) %>%
  dplyr::mutate(TCGA_barcode = gsub("\\.","-", TCGA_barcode),
                TCGA_barcode = TCGAbiospec(TCGA_barcode)$submitter_id)
cell_types <- colnames(rel)[-1]
combined_data <- merge(rel, survival_data, by.x = 'TCGA_barcode', by.y = 'bcr_patient_barcode', all.x = TRUE) %>%
  dplyr::filter(type %in% total_cancers) %>%
  dplyr::mutate(type = factor(type, levels = total_cancers))
combined_data <- split(combined_data, combined_data$type)
all(as.character(unname(unlist(lapply(combined_data, function(x) unique(x$type))))) == total_cancers)

# run analysis
for (i in seq_along(total_cancers)) {
  res <- pbmclapply(cell_types, cibersortlm, combined_data[[i]], conf = confounders, mc.cores = 10)
  if(i == 1){
    lm.res = lapply(res, pluck, 1) %>% bind_rows()
    cor.res = lapply(res, pluck, 2) %>% bind_rows()
  } else{
    lm.res = bind_rows(lm.res, lapply(res, pluck, 1) %>% bind_rows())
    cor.res = bind_rows(cor.res, lapply(res, pluck, 2) %>% bind_rows())
  }
}

cibersort.rel.lm <- lm.res %>%
  dplyr::filter(term == "QuartilesOld") %>%
  dplyr::group_by(cancer) %>%
  dplyr::mutate(fdr = p.adjust(p.value, "fdr"),
                dir = ifelse(estimate < 0 & fdr < 0.05, "Younger",
                             ifelse(estimate > 0 & fdr < 0.05, "Older",
                                    "NS")),
                dir = factor(dir, levels = c("Younger", "Older", "NS")),
                variable = paste0(cancer,"-",cell_type),
                cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  filter(fdr < 0.05)

cibersort.rel.corr <- cor.res %>%
  dplyr::group_by(cancer) %>%
  dplyr::mutate(fdr = p.adjust(p, "fdr"),
                dir = ifelse(cor < 0 & fdr < 0.05, "Younger",
                             ifelse(cor > 0 & fdr < 0.05, "Older",
                                    "NS")),
                dir = factor(dir, levels = c("Younger", "Older", "NS")),
                variable = paste0(cancer,"-",cell_type),
                cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated")) %>%
  dplyr::filter(fdr < 0.05)

## CIBERSORT RELATIVE + XCELL ##
################################
# linear model
bind_rows(cibersort.rel.lm %>% dplyr::mutate(tool  = "CIBERSORT Relative"),
          xcell.lm %>% dplyr::mutate(tool = "xCell")) %>%
  mutate(dir = factor(ifelse(dir == "Younger", "Young", "Old"), levels = c("Young", "Old"))) %>%
  filter(fdr < 0.05, cancer %in% naa_cancers) %>%
  ggplot(., aes(y = reorder(variable, estimate), x = estimate,
                xmin = estimate-std.error, xmax = estimate+std.error,
                color = dir, shape = cansource)) +
  facet_wrap(~tool, scales = "free_y") +
  geom_point(size=5) +
  geom_errorbar(size=2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Regression Coefficient", y = NULL, 
       color = "Association", caption = "FDR < 0.01",
       shape = "Cancer Type") +
  scale_color_jco() + xlim(-1.2,1.2) +
  theme(legend.position = 'top', legend.box="vertical")
ggsave("results/figures/fig2/xcell_and_cibersort_relative_fdr01.eps", width = 15)

# correlation
bind_rows(cibersort.rel.corr %>% dplyr::mutate(tool  = "CIBERSORT Relative"),
          xcell.corr %>% dplyr::mutate(tool = "xCell")) %>%
  mutate(dir = factor(ifelse(dir == "Younger", "Young", "Old"), levels = c("Young", "Old"))) %>%
  filter(fdr < 0.05) %>%
  ggplot(., aes(y = reorder(variable, cor), x = cor,
                color = dir)) +
  facet_wrap(~tool, scales = "free_y") +
  geom_point(size = 4) +
  labs(x = "Spearman", y = NULL, 
       caption = "FDR < 0.01") +
  scale_color_jco() + xlim(-0.35,0.35) + theme_bw(20) +
  theme(legend.position = 'none')
ggsave("results/figures/fig2/xcell_and_cibersort_relative_spearnab_fdr01.eps", width = 12, height = 10)


bind_rows(cibersort.rel.lm %>% dplyr::mutate(tool  = "CIBERSORT Relative"),
          xcell.lm %>% dplyr::mutate(tool = "xCell")) %>%
  mutate(dir = factor(ifelse(dir == "Younger", "Young", "Old"), levels = c("Young", "Old"))) %>%
  filter(fdr < 0.05, cancer %in% aa_cancers) %>%
  ggplot(., aes(y = reorder(variable, estimate), x = estimate,
                xmin = estimate-std.error, xmax = estimate+std.error,
                color = dir)) +
  facet_wrap(~tool, scales = "free_y") +
  geom_point(size=4) +
  geom_errorbar(size=2, width = 0.5) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Regression Coefficient", y = NULL, 
       color = "Association") +
  scale_color_jco() + xlim(-1.2,1.2) + theme_bw(20) +
  theme(legend.position = 'none', legend.box="vertical")
ggsave("results/figures/fig2/xcell_and_cibersort_relative_regression_aa_fdr05.eps", width = 15, height = 10)

bind_rows(cibersort.rel.lm %>% dplyr::mutate(tool  = "CIBERSORT Relative"),
          xcell.lm %>% dplyr::mutate(tool = "xCell")) %>%
  mutate(dir = factor(ifelse(dir == "Younger", "Young", "Old"), levels = c("Young", "Old"))) %>%
  filter(fdr < 0.05, cancer %in% naa_cancers) %>%
  ggplot(., aes(y = reorder(variable, estimate), x = estimate,
                xmin = estimate-std.error, xmax = estimate+std.error,
                color = dir)) +
  facet_wrap(~tool, scales = "free_y") +
  geom_point(size=5) +
  geom_errorbar(size=2) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Regression Coefficient", y = NULL) +
  scale_color_jco() + xlim(-1.2,1.2) + theme_bw(20) + 
  theme(legend.position = 'none', legend.box="vertical")
ggsave("results/figures/fig2/xcell_and_cibersort_relative_regression_naa_fdr05.eps", width = 7, height = 10)


## CIBERSORT ABSOLUTE ##
########################
cols <- setNames(pal_jco()(3), c("Young", "Old", "NS"))

# ## linear model
# # plot aa cancers
# cibersort.abs.lm %>%
#   mutate(cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
#          dir = ifelse(estimate < 0 & p.value.adj < 0.05, "Young",
#                       ifelse(estimate > 0 & p.value.adj < 0.05, "Old", "NS")),
#          dir = factor(dir, levels = c("Young", "Old", "NS"))) %>%
#   filter(cansource == "Age-Associated") %>%
#   ggplot(., aes(y = reorder(cancer, estimate), x = estimate,
#                 xmin = estimate-std.error, xmax = estimate+std.error,
#                 color = dir)) +
#   geom_point(size=5) +
#   geom_errorbar(size=2) +
#   geom_vline(xintercept = 0, lty = 2) +
#   labs(x = "Regression Coefficient", y = NULL, 
#        color = "Association", caption = "FDR < 0.05",
#        shape = "Cancer Type") +
#   scale_color_manual(values = cols) + xlim(-0.66, 0.66) +
#   theme(legend.position = 'top', legend.box="vertical")
# ggsave("results/figures/fig2/cibersort_abs_aa_fdr5.eps")
# 
# # plot naa cancers
# cibersort.abs.lm %>%
#   mutate(cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
#          dir = ifelse(estimate < 0 & p.value.adj < 0.05, "Young",
#                       ifelse(estimate > 0 & p.value.adj < 0.05, "Old", "NS")),
#          dir = factor(dir, levels = c("Young", "Old", "NS"))) %>%
#   filter(cansource == "Not Age-Associated") %>%
#   ggplot(., aes(y = reorder(cancer, estimate), x = estimate,
#                 xmin = estimate-std.error, xmax = estimate+std.error,
#                 color = dir)) +
#   geom_point(size=5) +
#   geom_errorbar(size=2) +
#   geom_vline(xintercept = 0, lty = 2) +
#   labs(x = "Regression Coefficient", y = NULL, 
#        color = "Association", caption = "FDR < 0.05",
#        shape = "Cancer Type") +
#   scale_color_manual(values = cols) + xlim(-0.66, 0.66) +
#   theme(legend.position = 'top', legend.box="vertical")
# ggsave("results/figures/fig2/cibersort_abs_aa_fdr5.eps")


## correlation
# cibersort.abs.corr %>%
#   mutate(cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
#          dir = ifelse(cor < 0 & p < 0.01, "Young",
#                       ifelse(cor > 0 & p < 0.01, "Old", "NS")),
#          dir = factor(dir, levels = c("Young", "Old", "NS"))) %>%
#   filter(cansource == "Age-Associated") %>%
#   ggplot(., aes(y = reorder(cancer, cor), x = cor, fill = dir)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Spearman", y = NULL, 
#        fill = NULL, caption = "FDR < 0.01",
#        shape = "Cancer Type") +
#   scale_fill_manual(values = cols)  + xlim(-0.2,0.2) +
#   theme(legend.position = 'top', legend.box="vertical")
# ggsave("results/figures/fig2/cibersort_abs_spearman_aa_fdr1.eps")
# 
# cibersort.abs.corr %>%
#   mutate(cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
#          dir = ifelse(cor < 0 & p < 0.01, "Young",
#                       ifelse(cor > 0 & p < 0.01, "Old", "NS")),
#          dir = factor(dir, levels = c("Young", "Old", "NS"))) %>%
#   filter(cansource == "Not Age-Associated") %>%
#   ggplot(., aes(y = reorder(cancer, cor), x = cor, fill = dir)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Spearman", y = NULL, 
#        fill = NULL, caption = "FDR < 0.01",
#        shape = "Cancer Type") +
#   scale_fill_manual(values = cols)  + xlim(-0.2,0.2) +
#   theme(legend.position = 'top', legend.box="vertical")
# ggsave("results/figures/fig2/cibersort_abs_spearman_naa_fdr1.eps")

cibersort.abs.corr %>%
  mutate(cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
         dir = ifelse(cor < 0 & p.adj < 0.05, "Young",
                      ifelse(cor > 0 & p.adj < 0.05, "Old", "NS")),
         dir = factor(dir, levels = c("Young", "Old", "NS"))) %>%
  # filter(cansource == "Age-Associated") %>%
  ggplot(., aes(y = reorder(cancer, cor), x = cor, color = dir, shape = cansource)) +
  geom_point(size = 4) +
  labs(x = "Spearman", y = NULL, 
       fill = NULL, caption = "FDR < 0.01",
       shape = "Cancer Type", color = "Correlation") +
  scale_color_manual(values = cols)  + xlim(-0.2,0.2) +
  theme_bw(20) +
  theme(legend.position = 'top', legend.box="vertical")
ggsave("results/figures/fig2/cibersort_abs_spearman_fdr5.eps", width = 7, height = 10)

cibersort.abs.lm %>%
  mutate(cansource = ifelse(cancer %in% aa_cancers, "Age-Associated", "Not Age-Associated"),
         dir = ifelse(estimate < 0 & p.value.adj < 0.05, "Young",
                      ifelse(estimate > 0 & p.value.adj < 0.05, "Old", "NS")),
         dir = factor(dir, levels = c("Young", "Old", "NS"))) %>%
  ggplot(., aes(y = reorder(cancer, estimate), x = estimate,
                color = dir,xmin = estimate-std.error, xmax = estimate+std.error)) +
  geom_point(size = 4) + geom_errorbar(size = 2, width = 0.5) + scale_color_manual(values = cols) + xlim(-0.75, 0.75) +
  labs(x = 'Regression Coefficient', y = NULL, color = 'Association', shape = 'Cancer Type') + facet_wrap(~cansource, scales = 'free_y', nrow = 2) +
  theme_bw(20) + geom_vline(xintercept = 0, lty=2) + theme(legend.position = 'none', legend.box = 'vertical')
ggsave("results/figures/fig2/cibersort_abs_regressopm_fdr5.eps", width = 12, height = 14)


## ICB signature
source("code/useful_functions.R")
source("code/02 - survival.R")
deg = readRDS("results/rds/archived/logistic/deg_gene.rds")[c(aa_cancers, naa_cancers)]
icb.res = mclapply(deg, function(this.df){
  dat = this.df %>%
    rownames_to_column('gene') %>%
    mutate(aa = ifelse(type %in% aa_cancers, 'Age-Associated', 'Not Age-Associated'))
  runGSEA(diff.results = dat, pathways = list(ICB = icb.genes)) %>% mutate(type = dat$type[1], aa = dat$aa[1])
}, mc.cores = 8)
rbindlist(icb.res) %>% filter(pval < 0.05) %>%
  ggplot(aes(x = NES, y = reorder(type, NES), fill = aa)) +
  geom_col() + scale_fill_manual(values = cols) +
  labs(x = 'NES', y = NULL, fill = NULL) +
  theme_pubr() +
  theme(legend.position = c(0.2, 0.9)) 
ggsave("results/figures/fig2/27_gene_icb_signature.eps", width = 6, height = 6.5)
