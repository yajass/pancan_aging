library(data.table)
library(tidyverse)
library(survival)
library(ggpubr)

survival_data <- openxlsx::read.xlsx("results/spreadsheets/survival_data.xlsx")
confounders <- openxlsx::read.xlsx("results/spreadsheets/confounder_identification/confounders_to_use.xlsx")
confounders[is.na(confounders)] = ""

# remove NA OS instances/time
vars <- data.frame(table(survival_data$type)) %>% filter(Freq > 100) %>% pull(Var1) %>% as.character(.)
vars <- vars[!vars %in% c("SKCM")] # less than 100 gene expression samples in these types
setdiff(survival_data$type %>% unique %>% as.character, vars)  ## removed due to lack of sample size
survival_data <- survival_data %>% filter(type %in% vars)
survival_data <- survival_data[!is.na(survival_data$OS),]
survival_data <- survival_data[!is.na(survival_data$OS.time),]
confounders <- confounders %>% filter(type %in% vars)
all(unique(survival_data$type) == confounders$type)

# loop to calculate OS survival
cancer_code <- unique(survival_data$type)
results_df <- data.frame(matrix(nrow = length(cancer_code),
                                ncol = 7))
colnames(results_df) <- c("Cancer", "Surv_N", "Surv_events", "hazard_ratio",
                          "lower_ci", "upper_ci", "Surv_pval")

for (i in 1:length(cancer_code)) {
  cancerDF <- survival_data[survival_data$type == cancer_code[i], ]
  
  # prep confounders for model matrix
  conf.vars <- confounders$linear[i]
  if (nchar(conf.vars) > 2){
    if(stringr::str_detect(conf.vars, pattern = "\\+")){
      conf.vars <- unlist(stringr::str_split(conf.vars, "\\+"))
    }
  }
  if (nchar(conf.vars[1]) < 2){
    conf.vars <- NULL
  }
  
  cancerDF <- cancerDF %>%
    dplyr::select(OS, OS.time, age_at_initial_pathologic_diagnosis, conf.vars)
  
  # coxph
  vars <- setdiff(colnames(cancerDF), c("OS","OS.time"))
  fS <- Surv(OS.time, OS) ~ . 
  fs <- reformulate(vars, fS[[2]])
  surv_result <- summary(coxph(fs, data = cancerDF))
  
  results_df$Surv_pval[i] <- coef(surv_result)[1,5]
  results_df$hazard_ratio[i] = round(surv_result$coefficients[1,2],5)
  results_df$Surv_N[i] <- surv_result$n
  results_df$Surv_events[i] <- surv_result$nevent
  results_df$lower_ci[i] <- surv_result$conf.int[1,3]
  results_df$upper_ci[i] <- surv_result$conf.int[1,4]
  results_df$Cancer[i] <- cancer_code[i]
}

# FDR
sig.level <- 0.05
results_df <- results_df %>%
  dplyr::mutate(fdr = p.adjust(results_df$Surv_pval, "fdr"),
                psig = ifelse(Surv_pval < 0.05, TRUE, FALSE),
                fdrsig = ifelse(fdr < sig.level, TRUE, FALSE),
                siganno = ifelse(psig & fdrsig, paste0("FDR < ", sig.level), "NS"),
                siganno = factor(siganno, levels = c(paste0("FDR < ", sig.level),"NS"))) %>%
  dplyr::arrange(siganno, -hazard_ratio)
levs <- rev(results_df$Cancer)
results_df$Cancer <- factor(results_df$Cancer, levels = levs)

results_df1 <- results_df %>% mutate(Cancer = paste0(Cancer, ' (n = ', Surv_N, ')'))  %>%
  arrange(fdrsig, hazard_ratio) %>% mutate(Cancer = factor(Cancer, levels = Cancer))

plt1 <- ggplot(data = results_df1, aes(y = Cancer, x = hazard_ratio, xmin = lower_ci, xmax = upper_ci, color = siganno)) +
  geom_pointrange(shape = 18, size=2) +
  geom_errorbar(width = 0.75, lwd = 1.5) +
  geom_vline(xintercept=1, lty=2, size = 2) +
  ggsci::scale_color_jama() +
  labs(x = "Hazard Ratio", y = NULL, title = NULL, color = NULL) +
  xlim(c(0.75, 1.25)) +
  theme_pubclean(20) +
  theme(axis.text.y = element_text(),
        legend.direction = "vertical",
        legend.justification = c("left", "top"),
        legend.position = c(.05, .95),
        axis.text.x = element_text())
plt1
# ggsave("results/figures/fig1/fig.1a.eps", width = 9, height = 8); rm(results_df1)

ggscatter(results_df, x = 'Surv_N', y = 'hazard_ratio', add = 'reg.line', color = 'siganno', palette = 'jama', ) +
  ggrepel::geom_text_repel(aes(label = Cancer, color = siganno), show.legend = F) + stat_cor(label.y.npc = 0.95, aes(color = siganno), show.legend = F) +
  stat_regline_equation(label.y.npc = 0.85, aes(color = siganno), show.legend = F) + theme_pubr(20) +
  labs(x = 'Sample Size', y = 'Harard Ratio', color = 'Hazard Ratio FDR')

levs = survival_data %>%
  # filter(Quartiles %in% c(1,4)) %>%
  # mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', 'Old'), levels = c('Young', 'Old'))) %>%
  group_by(type) %>%
  summarise(median = age_at_initial_pathologic_diagnosis %>% median)%>% arrange(-median) %>% pull(type)

survival_data %>% 
  filter(Quartiles %in% c(1,4)) %>%
  mutate(Quartiles = factor(ifelse(Quartiles == 1, 'Young', 'Old'), levels = c('Young', 'Old')),
         type = factor(type, levels = levs)) %>%
  ggplot(aes(x = Quartiles, y = age_at_initial_pathologic_diagnosis, fill = Quartiles)) + geom_violin(trim = F) +
  geom_boxplot(width = 0.2) + facet_wrap(~type, nrow = 1) + scale_fill_jco() +
  labs(y = 'Age') +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_text(size = 10)
  )
# ggsave('results/figures/fig1/quartile_age.eps', width = 18, height = 4)

sig_cancers <- results_df %>% filter(siganno == "FDR < 0.05") %>% pull(Cancer) %>% as.character(.)

survival_data1 <- survival_data %>% mutate(followup = ifelse(is.na(last_contact_days_to), death_days_to, last_contact_days_to),
                                           followup = as.numeric(as.character(followup))) %>%
  group_by(type)

survival_data1 = survival_data1 %>%
  rstatix::cor_test(vars = age_at_initial_pathologic_diagnosis, vars2 = followup, method = 'pearson') %>%
  rstatix::adjust_pvalue(method = 'fdr') %>%
  mutate(col = ifelse(p.adj < 0.05, "FDR < 0.05", "NS"))

survival_data1 %>%
  ggplot(aes(x=cor,y=-log10(p.adj), label=type, color = col)) +
  geom_point() + 
  ggrepel::geom_label_repel(data = subset(survival_data1, col != 'NS'), show.legend = F) + 
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  ggsci::scale_color_jama() + labs(x = "Pearson", y = "-log FDR", color = NULL) +
  theme_bw(20) + theme(legend.position = 'top')
## ggsave('results/figures/fig1/followup_time_correlatoion.eps', width = 5, height = 5)

rm(survival_data1)

# fwrite(results_df, 'results/tables/coxph_output.tsv', row.names = F, sep = '\t', quote = F)

# survival_data %>%
#   filter(Quartiles %in% c(1, 4)) %>%
#   mutate(Age = ifelse(Quartiles == 1, 'Young', 'Old')) %>%
#   group_by(type, Age) %>%
#   summarize(Minimum = min(age_at_initial_pathologic_diagnosis),
#             Maximum = max(age_at_initial_pathologic_diagnosis),
#             Mean = mean(age_at_initial_pathologic_diagnosis),
#             Median = median(age_at_initial_pathologic_diagnosis)) %>%
#   fwrite('results/tables/age_classification.tsv', row.names = F, sep = '\t', quote = F)
