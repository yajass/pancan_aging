library(data.table)
library(tidyverse)
library(survival)
library(ggpubr)

survival_data <- openxlsx::read.xlsx("results/spreadsheets/survival_data.xlsx")
confounders <- openxlsx::read.xlsx("results/spreadsheets/confounder_identification/confounders_to_use.xlsx")
confounders[is.na(confounders)] = ""

# remove NA OS instances/time
vars <- data.frame(table(survival_data$type)) %>% filter(Freq > 100) %>% pull(Var1) %>% as.character(.)
vars <- vars[!vars %in% c("SKCM")] # less than 150 gene expression samples in these types
survival_data <- survival_data %>% filter(type %in% vars)
survival_data <- survival_data[!is.na(survival_data$OS),]
survival_data <- survival_data[!is.na(survival_data$OS.time),]
confounders <- confounders %>% filter(type %in% vars)
all(unique(survival_data$type) == confounders$type)

cancer_code <- unique(survival_data$type)


doSurvival <- function(cancer_code, survival_data, confounders){
  cancerDF1 <- survival_data[survival_data$type == cancer_code, ]
  
  idx <- sample(1:nrow(cancerDF1), 100, replace = FALSE)
  cancerDF <- cancerDF1[idx, ]
  
  while (length(table(cancerDF$OS)) != 2) {
    idx <- sample(1:nrow(cancerDF1), 100, replace = FALSE)
    cancerDF <- cancerDF1[idx, ]
  }
  
  # prep confounders for model matrix
  conf.vars <- confounders$linear[confounders$type == cancer_code]
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
  
  keep <- names(apply(cancerDF, 2, function(x) table(x) %>% length) > 1)
  
  cancerDF <- cancerDF %>%
    dplyr::select(all_of(keep))
  
  # coxph
  vars <- setdiff(colnames(cancerDF), c("OS","OS.time"))
  fS <- Surv(OS.time, OS) ~ . 
  fs <- reformulate(vars, fS[[2]])
  surv_result <- summary(coxph(fs, data = cancerDF))
  
  # reutrn results
  results_df <- data.frame(matrix(nrow = 1,
                                  ncol = 7))
  colnames(results_df) <- c("Cancer", "Surv_N", "Surv_events", "hazard_ratio",
                            "lower_ci", "upper_ci", "Surv_pval")
  results_df$Surv_pval <- coef(surv_result)[1,5]
  results_df$hazard_ratio = round(surv_result$coefficients[1,2],5)
  results_df$Surv_N <- surv_result$n
  results_df$Surv_events <- surv_result$nevent
  results_df$lower_ci <- surv_result$conf.int[1,3]
  results_df$upper_ci <- surv_result$conf.int[1,4]
  results_df$Cancer <- cancer_code
  
  return(results_df)
}

for (rep in 1:25) {
  if (rep == 1){
    res = pbmclapply(cancer_code, doSurvival, survival_data = survival_data, confounders = confounders, mc.cores=13) %>% bind_rows() %>% mutate(rep = rep)
  }
  if (rep > 1){
    res = bind_rows(res,
                    pbmclapply(cancer_code, doSurvival, survival_data = survival_data, confounders = confounders, mc.cores=13) %>% bind_rows() %>% mutate(rep = rep))
  }
}

res1 <- res %>%
  group_by(rep) %>%
  mutate(fdr = p.adjust(Surv_pval, 'fdr'),
         sig = factor(ifelse(fdr < 0.05, "FDR < 0.05", "NS"), levels = c("NS", "FDR < 0.05"))) %>% 
  ungroup()

dplyr::count(res1, Cancer, sig) %>% arrange(sig, -n) %>%
  mutate(Cancer = factor(Cancer, levels = unique(Cancer))) %>%
  ggbarplot(x = 'Cancer', y = 'n', fill = 'sig', palette = 'Set1') +
  labs(x = NULL, y = 'No. of Replicates', fill=  NULL) + theme_minimal(20) +
  coord_flip() + theme(legend.position = 'top')
ggsave('results/figures/fig1/survival_validation.eps', width = 5, height = 7)
dev.off()

tumors = dplyr::count(res1, Cancer, sig) %>% arrange(sig, -n) %>%
  mutate(Cancer = factor(Cancer, levels = unique(Cancer))) %>%
  filter(sig == 'FDR < 0.05') %>% pull(Cancer) %>% as.character

intersect(sig_cancers, tumors)

cairo_ps('results/figures/fig1/survival_validation_venn.eps')
ggvenn::ggvenn(data = list(`Without\nDownsampling` = sig_cancers,
                           `With\nDownsampling` = tumors),
               show_elements = T, label_sep = "\n", stroke_linetype = 0, text_size = 2.5)
dev.off()

saveRDS(res, 'results/rds/downsampling_survival.rds')
