library(ggsci); library(ggpubr); theme_set(theme_pubr()); options(stringsAsFactors = FALSE)
library(tidyverse); library(data.table); library(tidyverse)
library(TCGAbiolinks); library(TCGAutils)
source("code/useful_functions.R")
source("code/02 - survival.R")

sv <- survival_data %>% dplyr::select(bcr_patient_barcode, type, Quartiles)
dir <- "../Data/methylation_450/horvath_clock/output/results/"

clock <- lapply(list.files(dir), function(x) fread(paste0(dir, x))) %>% 
  bind_rows() %>%
  dplyr::mutate(definition = TCGAbiospec(gsub("\\.","-",SampleID))$sample_definition) %>%
  inner_join(sv) %>%
  filter(corSampleVSgoldstandard >= 0.8,
         definition %in% c("Primary Solid Tumor", "Primary Blood Derived Cancer - Peripheral Blood", "Solid Tissue Normal"),
         type %in% c(aa_cancers, naa_cancers)) %>%
  mutate(definition = ifelse(definition == "Solid Tissue Normal", "Normal", "Tumor"),
         definition = factor(definition, levels = c("Tumor", "Normal")),
         class = ifelse(type %in% aa_cancers, "AA", "NAA")) %>%
  arrange(rev(definition))
paired_samples <- clock %>% count(bcr_patient_barcode) %>% filter(n > 1) %>% pull(bcr_patient_barcode)
clock %>%
  dplyr::select(SampleID, type, definition, Age, AgeAccelerationDiff, DNAmAge, 
                corSampleVSgoldstandard, meanAbsDifferenceSampleVSgoldstandard, class) %>%
  fwrite('results/tables/epigenetic_age.tsv', row.names = F, sep = '\t', quote = F)


# look for outliers
clock %>% filter(type != 'OV') %>%
  ggscatter(x = 'Age', y = 'DNAmAge', facet.by = c('definition', 'type'), alpha = 2/3) + stat_cor(label.sep = '\n') +
  theme_bw(25)
ggsave("results/figures/fig3/dnamage_qc.jpg", width = 28, height = 8)

# Plot age associated
plt3a = ggscatter(clock %>% filter(class == "AA", type != "OV"), x = "Age", y = "DNAmAge", color = "definition", add = "reg.line",
          conf.int = TRUE, alpha = 1/3) + labs(color = NULL, y = "DNAm Age") + scale_color_manual(values = cols) +
  stat_cor(aes(color = definition), show.legend = FALSE) + theme_pubr(25) +
  stat_regline_equation(aes(color = definition), label.y.npc = 0.9, show.legend = F) +
  guides(fill = FALSE) + theme(legend.position = c(0.55, 0.95), legend.direction = 'horizontal')
## ggsave("results/figures/fig3/unpaired_aa_dnamage.jpg", width = 8, height = 7)

# Plot NAA
ggscatter(clock %>% filter(class == "NAA", type != "OV"), x = "Age", y = "DNAmAge", color = "definition", add = "reg.line",
          conf.int = TRUE, alpha = 1/3) + labs(color = NULL, y = "DNAm Age") + scale_color_manual(values = cols) +
  stat_cor(aes(color = definition), show.legend = FALSE) + theme_pubr(20) +
  stat_regline_equation(aes(color = definition), label.y.npc = 0.9, show.legend = F) +
  guides(fill = FALSE) + theme(legend.position = "top")
ggsave("results/figures/fig3/unpaired_naa_dnamage.jpg", width = 8, height = 7)

# Age acceleration difference - unpaired
clk = clock %>% filter(definition == 'Tumor', type != 'OV') %>% filter(Quartiles %in% c(1,4)) %>%  
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', 'Old'), levels = c('Young', 'Old')))
st = clk %>% group_by(Quartiles) %>%
  rstatix::wilcox_test(AgeAccelerationDiff  ~ class) %>% rstatix::add_significance()
ggviolin(clk, x = 'class', y = 'AgeAccelerationDiff', facet.by = 'Quartiles', add = 'boxplot', fill = 'class') + theme_pubr(20) +
  scale_fill_manual(values = cols) + theme(legend.position = 'none') + labs(x = NULL, y = 'DNAm Age Acc Diff') +
  stat_pvalue_manual(data = st, y = 75, bracket.size = 0, size = 10)
ggsave("results/figures/fig3/dnam_tumor_aa_naa.jpg", width = 7, height = 3.5)

# plot decreased age acc diff with aging in tumor and normal
st = clock %>% filter(type != 'OV') %>%
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', ifelse(Quartiles == 4, 'Old', 'Middle\nAged')), levels = c('Young', 'Middle\nAged', 'Old'))) %>% 
  group_by(definition) %>%
  rstatix::wilcox_test(AgeAccelerationDiff ~ Quartiles) %>% filter(group1 == 'Young', group2 == 'Old')
clock %>% filter(type != 'OV') %>%
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', ifelse(Quartiles == 4, 'Old', 'Middle\nAged')), levels = c('Young', 'Middle\nAged', 'Old'))) %>%
  ggviolin(x = 'Quartiles', y = 'AgeAccelerationDiff', fill = 'class', add = 'boxplot', facet.by = 'definition') +
  labs(x = NULL, y = 'DNAm Age Acc Diff', fill = NULL) + expand_limits(y = 130) + theme_pubr(22, legend = 'top') + scale_fill_manual(values = cols) +
  stat_compare_means(aes(group = class), method = 'wilcox', label.y = 95, label = 'p.signif', size = 10) +
  stat_pvalue_manual(data = st, y.position = 120, size = 10, bracket.size = 1.2)
ggsave('results/figures/fig3/decreased_age_acc_violin.jpg', width = 8, height = 7)


plt3b = clock %>% filter(type != 'OV', definition == 'Tumor') %>% arrange(class) %>%
  ggscatter(x = 'Age', y = 'AgeAccelerationDiff', color = 'class', alpha = 1/3) + geom_smooth(aes(color = class), method = 'lm', se = T) +
  scale_color_manual(values = cols) + stat_cor(aes(color = class), show.legend = F) +
  stat_regline_equation(aes(color = class), show.legend = F, label.y.npc = 0.9) + theme_pubr(25) +
  labs(y = 'DNAm Age Acc Diff', color = NULL) + theme(legend.position = c(0.55,0.95), legend.direction = 'horizontal')
# ggsave('results/figures/fig3/dnam_age_acc_scatter.jpg', width = 9, height = 9, dpi = 320)

ggarrange(plt3a, plt3b, align = 'h')
ggsave('results/figures/fig3/dnam_3a_3b.jpg', width = 16, height = 7, dpi = 320)




# Plot age acceleration diff delta (tumor - normal) acros age assoc and non age assoc
clock %>%
  dplyr::select(bcr_patient_barcode, type, definition, class, AgeAccelerationDiff, Age) %>%
  filter(bcr_patient_barcode %in% paired_samples, type != "OV") %>%
  group_by(type) %>% mutate(scaled = (AgeAccelerationDiff - mean(AgeAccelerationDiff))/sd(AgeAccelerationDiff)) %>% ungroup %>%
  group_by(type, bcr_patient_barcode, class, Age) %>%
  summarise(delta_t_n = scaled[definition == "Tumor"] - scaled[definition == "Normal"]) %>%
  ggscatter(., x = "Age", y = "delta_t_n", color = "class", add = "reg.line", conf.int = TRUE, alpha = 3/4) + scale_color_manual(values = cols) +
  stat_cor(aes(color = class), show.legend = F) + stat_regline_equation(aes(color = class), show.legend = F, label.y.npc = 0.9) +
  theme_pubr(20) + labs(y = "∆ Scaled AgeAccDiff (Tumor-Normal)", color = NULL) + 
  guides(fill = FALSE) + theme(legend.position = "top")
ggsave("results/figures/fig3/paired_aa_vs_naa_dnamage.jpg")

# Do statistics on age acceleration diff delta (tumor - normal) acros age assoc and non age assoc
p <- clock %>%
  dplyr::select(bcr_patient_barcode, type, definition, class, AgeAccelerationDiff, Age, Quartiles) %>%
  group_by(Quartiles, type, bcr_patient_barcode, class, Age) %>%
  filter(bcr_patient_barcode %in% paired_samples, Quartiles %in% c(1,4), type != "OV") %>%
  summarise(delta_t_n = AgeAccelerationDiff[definition == "Tumor"] - AgeAccelerationDiff[definition == "Normal"]) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"),
         Quartiles = factor(Quartiles, levels = c("Young", "Old"))) %>%
  ggviolin(., x = "Quartiles", y = "delta_t_n", facet.by = "class", fill = "Quartiles", palette = "jco", add = "boxplot") +
  stat_compare_means(method = "t.test", comparisons = list(c("Young","Old")), label.y = 40, label = "p.signif", size = 5, bracket.size = 0, hide.ns = T) +
  labs(y = "∆ Age Acc Diff (T-N)", fill = "Age", x = NULL) + theme_pubr(20, legend = 'none')
p$layers[[3]]$aes_params$textsize <- 10; p
ggsave("results/figures/fig3/paired_aa_vs_naa_dnamage_violin.jpg", width = 7, height = 3.5)

# plot tumor types contibuting to age acceleration phenotype - ONLY TUMOR SAMPLE
p <- clock %>%
  dplyr::filter(class == "AA", definition == "Tumor", type != "OV") %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", "Middle Aged")),
         Quartiles = factor(Quartiles, levels = c("Young", "Middle Aged", "Old"))) %>%
  ggplot(., aes(x = Quartiles, y = AgeAccelerationDiff, color = type, group = type)) +
  geom_errorbar(stat = 'summary', fun.data = 'mean_ci', width = 0.2, size = 1) +
  geom_point(stat = 'summary', fun = 'mean', size = 3) +
  stat_summary(fun = "mean", geom = "line", lty = 2, size = 1) +
  theme_pubr(20) + scale_color_npg() + expand_limits(y = 45) +
  labs(x = NULL, y = "DNAm AgeAccDiff", color = NULL)
  # stat_compare_means(comparisons = list(c("Young", "Old")), label.y = 40, methodd = "wilcox", label = "p.signif", bracket.size = 1.5, size = 10)
# p$layers[[2]]$aes_params$textsize <- 10; p
p
ggsave("results/figures/fig3/dnam_age_acc_by_type_aa.eps", width = 7, height = 3.5)


p <- clock %>%
  dplyr::filter(class == "NAA", definition == "Tumor", type != "OV") %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", ifelse(Quartiles == 4, "Old", "Middle Aged")),
         Quartiles = factor(Quartiles, levels = c("Young", "Middle Aged", "Old"))) %>%
  ggplot(., aes(x = Quartiles, y = AgeAccelerationDiff, color = type, group = type)) +
  geom_errorbar(stat = 'summary', fun.data = 'mean_ci', width = 0.2, size = 2) +
  geom_point(stat = 'summary', fun = 'mean', size = 4) +
  stat_summary(fun = "mean", geom = "line", lty = 2, size = 2) +
  theme_pubr(20) + scale_color_npg() + expand_limits(y = 35) +
  labs(x = NULL, y = "DNAm AgeAccDiff", color = NULL)
  # stat_compare_means(comparisons = list(c("Young", "Old")), label.y = 40, methodd = "wilcox", label = "p.signif", bracket.size = 1.5, size = 10)
# p$layers[[2]]$aes_params$textsize <- 10; p
p
ggsave("results/figures/fig3/dnam_age_acc_by_type_naa.eps", width = 8, height = 7)
