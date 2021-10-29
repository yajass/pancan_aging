library(maftools)
library(ggpubr); library(ggsci); library(gginnards)

# load data
source("code/useful_functions.R")
source("code/02 - survival.R")
aa_maf <- readRDS("results/rds/try.aa_maf.rds")
naa_maf <- readRDS("results/rds/naa_maf.rds")

# aa_maf = subsetMaf(aa_maf, clinQuery = "Dataset == 'TCGA'")
# plot TMB
dat <- aa_maf@clinical.data %>%
  dplyr::select(Quartiles, Tumor_Sample_Barcode, type, Dataset) %>%
  inner_join(aa_maf@variant.type.summary %>% dplyr::select(Tumor_Sample_Barcode, total)) %>%
  mutate(type = factor(type, levels = c("LUSC","UCEC","OV","BRCA","LGG","THCA")),
         Quartiles = factor(Quartiles, levels = c("Young","Old")),
         Dataset = factor(Dataset, levels = c('TCGA', 'METABRIC')))
stats <- dat %>%
  group_by(type, Dataset) %>%
  rstatix::wilcox_test(total ~ Quartiles) %>%
  rstatix::adjust_pvalue("p", "fdr", output.col = "p") %>% rstatix::add_significance(p.col = "p")
dat %>% mutate(total = total/3000) %>%
  ggviolin(., x = "Quartiles", y = "total", fill = "Quartiles", palette = "jco", add = 'boxplot') +
  facet_wrap(~ Dataset + type, nrow = 1) +
  scale_y_log10(name = 'TMB', breaks = c(0.01, 0.1, 1, 10, 100), labels = (c(0.01, 0.1, 1, 10, 100))) +
  stat_pvalue_manual(data = stats, label = "p.signif", size = 10, y.position = 1, hide.ns = T) +
  labs(x = NULL, fill = 'Age') + theme_pubr(20) + 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "lines")
  )
ggsave("results/figures/fig4/tmb_aa.eps", height = 7, width = 15)

# Non AA
## naa_maf = subsetMaf(naa_maf, clinQuery = "Dataset == 'TCGA'")
# plot TMB
dat <- naa_maf@clinical.data %>%
  dplyr::select(Quartiles, Tumor_Sample_Barcode, type) %>%
  inner_join(naa_maf@variant.type.summary %>% dplyr::select(Tumor_Sample_Barcode, total)) %>%
  mutate(Quartiles = factor(Quartiles, levels = c("Young","Old")),
         total = total/3000,
         type = as.character(type)) %>% ungroup()
levs = setDT(dat)[,list(median = median(total)), by=type] %>% arrange(-median) %>% pull(type)
dat$type = factor(dat$type, levels = levs)
stats <- dat %>%
  group_by(type) %>%
  rstatix::wilcox_test(total ~ Quartiles) %>%
  rstatix::adjust_pvalue("p", "fdr", output.col = "p") %>% rstatix::add_significance(p.col = "p")
dat %>%
  ggviolin(., x = "Quartiles", y = "total", fill = "Quartiles", palette = "jco", add = 'boxplot') +
  facet_wrap(~type, nrow = 1) +
  scale_y_log10(name = 'TMB', breaks = c(0.01, 0.1, 1, 10, 100), labels = (c(0.01, 0.1, 1, 10, 100))) +
  stat_pvalue_manual(data = stats, label = "p.signif", size = 10, y.position = 1, hide.ns = T) +
  labs(x = NULL, fill = 'Age') + theme_pubr(20) + 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "lines")
  )
ggsave("results/figures/fig4/tmb_naa.eps", height = 7, width = 17)

# plot pancan snv - AA
dat = aa_maf@clinical.data %>%
  dplyr::select(Quartiles, Tumor_Sample_Barcode, type) %>%
  inner_join(aa_maf@variant.classification.summary) %>%
  dplyr::select(-Tumor_Sample_Barcode, -total) %>%
  melt(c("Quartiles", "type")) %>%
  mutate(variable = gsub("\\_", " ", variable),
         Quartiles = factor(Quartiles, levels = c("Young","Old")))
stats <- dat %>%
  group_by(variable) %>%
  rstatix::t_test(value ~ Quartiles) %>%
  rstatix::adjust_pvalue("p", "fdr", output.col = "p") %>% rstatix::add_significance(p.col = "p")
dat %>%
  ggplot(., aes(x = Quartiles, y= value)) +
  facet_wrap(~variable, scales = "free_y") +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Quartiles)) + ylim(0, NA) +
  stat_summary(fun.data = "mean_ci", geom = "errorbar", width = 0.2, size = 2) +
  labs(x = NULL, y = "Mean frequency per patient") + theme_pubr(20, legend = 'none') + scale_fill_jco() +
  stat_pvalue_manual(stats, label = "p.signif",hide.ns = T, label.size = 10, bracket.size = 0, inheric_aes = T,
                     y.pos = c(280,28,0.3,7,0.28), tip.length = 0.0001)
ggsave("results/figures/fig4/aa_pancan_snvs.jpeg", width = 11, height = 6, dpi = 320)

# pancan SNV - NAA
dat = naa_maf@clinical.data %>%
  dplyr::select(Quartiles, Tumor_Sample_Barcode, type) %>%
  inner_join(naa_maf@variant.classification.summary) %>%
  dplyr::select(-Tumor_Sample_Barcode, -total) %>%
  melt(c("Quartiles", "type")) %>%
  mutate(variable = gsub("\\_", " ", variable),
         Quartiles = factor(Quartiles, levels = c("Young","Old")))
stats <- dat %>%
  group_by(variable) %>%
  rstatix::t_test(value ~ Quartiles) %>%
  rstatix::adjust_pvalue("p", "fdr", output.col = "p") %>% rstatix::add_significance(p.col = "p")
dat %>%
  ggplot(., aes(x = Quartiles, y= value)) +
  facet_wrap(~variable, scales = "free_y") +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Quartiles)) + ylim(0, NA) +
  stat_summary(fun.data = "mean_ci", geom = "errorbar", width = 0.2, size = 2) +
  labs(x = NULL, y = "Mean frequency per patient") + theme_pubr(20, legend = 'none') + scale_fill_jco() +
  stat_pvalue_manual(stats, label = "p.signif",hide.ns = T, label.size = 10, bracket.size = 0, inheric_aes = T,
                     y.pos = c(18, 4.5, 1.8, 0.25), tip.length = 0.0001)
ggsave("results/figures/fig4/naa_pancan_snvs.jpeg", width = 14, height = 5, dpi = 320)

# plot SNV types
dat <- aa_maf@clinical.data %>%
  dplyr::select(Quartiles, Tumor_Sample_Barcode, type) %>%
  inner_join(aa_maf@variant.classification.summary) %>%
  dplyr::select(-Tumor_Sample_Barcode, -total) %>%
  melt(c("Quartiles", "type")) %>%
  mutate(variable = gsub("\\_", " ", variable),
         Quartiles = factor(Quartiles, levels = c("Young","Old"))) %>%
  filter(type == "UCEC")
stats <- dat %>%
  group_by(variable) %>%
  rstatix::t_test(value ~ Quartiles) %>%
  rstatix::adjust_pvalue("p", "fdr", output.col = "p") %>% rstatix::add_significance(p.col = "p") %>% mutate(y.position = c(0, 0, 0, 0,2000,200, 2,45,1.6))
dat %>%
  ggplot(., aes(x = Quartiles, y= value)) +
  facet_wrap(~variable, scales = "free_y") +
  geom_bar(stat = "summary", fun = "mean", aes(fill = Quartiles)) +
  stat_summary(fun.data = "mean_ci", geom = "errorbar", width = 0.2, size = 2) +
  labs(x = NULL, y = "Mean frequency per patient") + theme_pubr(20, legend = 'none') + scale_fill_jco() +
  stat_pvalue_manual(stats, label = "p.signif", hide.ns = T, label.size = 10, bracket.size = 0, inheric_aes = T)
ggsave("results/figures/fig4/ucec_snvs.eps", width = 9, height = 9)
