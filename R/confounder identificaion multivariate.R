library(TCGAbiolinks)
library(TCGAutils)
library(tidyverse)
library(tidymodels)
library(data.table)
library(ggpubr)

# load data
followup <- fread("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Survival_Data/clinical_PANCAN_patient_with_followup.tsv") %>%
  dplyr::select(bcr_patient_barcode, tobacco_smoking_history, alcoholic_exposure_category) %>%
  dplyr::rename(smoking_status = tobacco_smoking_history,
                alcohol_status = alcoholic_exposure_category)
survival_data <- openxlsx::read.xlsx("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/Survival_Data/TCGA-CDR-SupplementalTableS1.xltx")
survival_data <- survival_data[!is.na(survival_data$age_at_initial_pathologic_diagnosis), ] %>%
  dplyr::mutate(type = ifelse(type == "READ", "COAD", type)) %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(Quartiles = dplyr::ntile(age_at_initial_pathologic_diagnosis, 4)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ajcc_pathologic_tumor_stage = ifelse(ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB", "IS"),
                                                     "Stage I",
                                                     ifelse(ajcc_pathologic_tumor_stage %in% c("Stage IIA", "Stage IIB", "Stage IIC"),
                                                            "Stage II",
                                                            ifelse(ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIB", "Stage IIC"),
                                                                   "Stage III",
                                                                   ifelse(ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"),
                                                                          "Stage IV",
                                                                          ifelse(clinical_stage %in% c("Stage IS", "I","Stage I","Stage IA","Stage IA1","Stage IA2","Stage IB","Stage IB1","Stage IB2","Stage IC"),
                                                                                 "Stage I",
                                                                                 ifelse(clinical_stage %in% c("IIa","IIb","Stage II","Stage IIA","Stage IIA1","Stage IIA2","Stage IIB","Stage IIC"),
                                                                                        "Stage II",
                                                                                        ifelse(clinical_stage %in% c("III","Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IIIC1","Stage IIIC2"),
                                                                                               "Stage III",
                                                                                               ifelse(clinical_stage %in% c("IVa","IVb","Stage IV","Stage IVA","Stage IVB","Stage IVC"),
                                                                                                      "Stage IV", "missing")))))))),
                histological_type = ifelse(is.na(histological_type), "missing", histological_type))
survival_data <- merge(survival_data, followup, all.x = TRUE) %>%
  dplyr::mutate(smoking_status = ifelse(smoking_status %in% c("", "[Discrepency]", "[Not Available]", "[Unknown]"), "[Not Available]",
                                        ifelse(smoking_status %in% c("Current reformed smoker for < or = 15 years", "Current reformed smoker for > 15 years","Current Reformed Smoker, Duration Not Specified"),
                                               "Reformed Smoked", smoking_status)),
                alcohol_status = ifelse(alcohol_status %in% c("","[Not Available]", "[Not Evaluated]", "[Unknown]"), "[Not Avalable]", alcohol_status))

survival_data_log <- survival_data %>%
  dplyr::arrange(type, bcr_patient_barcode) %>%
  dplyr::filter(Quartiles %in% c(1,4)) %>%
  dplyr::mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))
survival_data_lin <- survival_data %>%
  dplyr::arrange(type, bcr_patient_barcode)

# identify confounders using multivariate linear regression
cancers <- unique(survival_data_lin$type)
mod.summaries <- list()
colnames(survival_data_lin) <- paste0(colnames(survival_data_lin),"QQ")
colnames(survival_data_log) <- paste0(colnames(survival_data_log),"QQ")

for (i in 1:length(cancers)) {
  
  cancer_data <- survival_data_lin[,c(1,3,4,5,6,7,9,10,36,37)] %>%
    dplyr::filter(typeQQ == cancers[i]) %>%
    dplyr::select(-bcr_patient_barcodeQQ, -typeQQ)
  cancer_data <- data.frame(cancer_data[, apply(cancer_data, 2, function(x) length(table(x))) > 1])
    
  # run regression
  lm_fit <- lm(age_at_initial_pathologic_diagnosisQQ ~ ., data = cancer_data)
  
  if (i == 1){
    mod_results <- tidy(lm_fit) %>% dplyr::mutate(typeQQ = cancers[i])
  }
  if (i > 1){
    mod_results <- dplyr::bind_rows(mod_results,
                                    tidy(lm_fit) %>% dplyr::mutate(typeQQ = cancers[i]))
  }
  
}

mod_results1 <- mod_results
mod_results1$term <- stringr::str_replace(mod_results1$term, pattern = "QQ", replacement = "")
openxlsx::write.xlsx(mod_results1, "results/spreadsheets/confounder_identification/multivariate_linear_regression/model_summaries.xlsx",
                     asTable = T)
lin_mod <- mod_results %>% 
  dplyr::mutate(term = gsub( "QQ.*$", "", term))

# identify confounders using multivariate logistic regression
cancers <- unique(survival_data_log$typeQQ)
mod.summaries <- list()
for (i in 1:length(cancers)) {
  
  cancer_data <- survival_data_log[,c(1,3,35,5,6,7,9,10,36,37)] %>%
    dplyr::filter(typeQQ == cancers[i]) %>%
    dplyr::select(-bcr_patient_barcodeQQ, -typeQQ)
  cancer_data <- data.frame(cancer_data[, apply(cancer_data, 2, function(x) length(table(x))) > 1]) %>%
    dplyr::mutate(QuartilesQQ = factor(QuartilesQQ, levels = c("Young", "Old")))
  
  # run regression
  lm_fit <- arm::bayesglm(QuartilesQQ ~ ., family = "binomial", data = cancer_data)
  
  if (i == 1){
    mod_results <- tidy(lm_fit) %>% dplyr::mutate(typeQQ = cancers[i])
  }
  if (i > 1){
    mod_results <- dplyr::bind_rows(mod_results,
                                    tidy(lm_fit) %>% dplyr::mutate(typeQQ = cancers[i]))
  }
  
}

mod_results1 <- mod_results
mod_results1$term <- stringr::str_replace(mod_results1$term, pattern = "QQ", replacement = "")

openxlsx::write.xlsx(mod_results1, "results/spreadsheets/confounder_identification/multivariate_logistic_regression/model_summaries.xlsx",
                     asTable = T)

log_mod <- mod_results %>% 
  dplyr::mutate(term = gsub( "QQ.*$", "", term))


## Plot models
# plot linear
lin_mod <- lin_mod %>%
  dplyr::filter(term != "(Intercept)",
                p.value < 0.05) %>%
  dplyr::select(term, typeQQ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(value = "linear")
lin_mod %>%
  dplyr::mutate(value = "p < 0.05") %>%
  ggplot(., aes(x = typeQQ, y = term, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = ggsci::pal_aaas()(3)[1]) +
  theme_pubr() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/spreadsheets/confounder_identification/multivariate_linear_regression/sig_features.eps")

# plot logistic
log_mod <- log_mod %>%
  dplyr::filter(term != "(Intercept)",
                p.value < 0.05) %>%
  dplyr::select(term, typeQQ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(value = "logistic")
log_mod %>%
  dplyr::mutate(value = "p < 0.05") %>%
  ggplot(., aes(x = typeQQ, y = term, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = ggsci::pal_aaas()(3)[2]) +
  theme_pubr() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/spreadsheets/confounder_identification/multivariate_logistic_regression/sig_features.eps")

# plot combined
combined <- rbind(lin_mod, log_mod) %>%
  dplyr::mutate(comb = paste(term, typeQQ),
                value = ifelse(value == "linear", "Continuous Age", "Quartiled Age"))
combined[combined$comb %in% (combined %>%
           group_by(comb) %>%
           dplyr::count() %>% 
           filter(n > 1) %>%
           dplyr::pull(comb)),
         3] <- "Continuous & Quartiled Age"
combined %>%
  dplyr::distinct() %>%
  dplyr::mutate(value = factor(value, levels = c("Quartiled Age", "Continuous Age", "Continuous & Quartiled Age"))) %>%
  ggplot(., aes(x = typeQQ, y = term, fill = value)) +
  geom_tile() +
  ggsci::scale_fill_aaas() +
  theme_pubr() +
  labs(x = NULL, y = NULL, fill = "Confounder") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/fig1/confounder_sig_features_linear_and_logistic.eps", width = 8, height = 6)

# save significant reults as xlsx
lin_mod %>%
  dplyr::group_by(typeQQ) %>%
  dplyr::summarise(linear = stringr::str_c(term, collapse = "+"))
log_mod %>%
  dplyr::group_by(typeQQ) %>%
  dplyr::summarise(logistic = stringr::str_c(term, collapse = "+"))

results <- merge(data.frame(type = cancers),
                 lin_mod %>%
                   dplyr::group_by(typeQQ) %>%
                   dplyr::summarise(linear = stringr::str_c(term, collapse = "+")) %>%
                   dplyr::rename(type = typeQQ) %>%
                   dplyr::ungroup(), all.x = TRUE) %>%
  dplyr::mutate(linear = ifelse(is.na(linear), "", linear))
results <- merge(results, log_mod %>%
                   dplyr::group_by(typeQQ) %>%
                   dplyr::summarise(logistic = stringr::str_c(term, collapse = "+")) %>%
                   dplyr::rename(type = typeQQ) %>%
                   dplyr::ungroup(), all.x = TRUE)  %>%
  dplyr::mutate(logistic = ifelse(is.na(logistic), "", linear))

openxlsx::write.xlsx(results, "results/spreadsheets/confounder_identification/confounders_to_use.xlsx")

# write survival data
colnames(survival_data_lin) <- gsub('.{2}$', '', colnames(survival_data_lin))
openxlsx::write.xlsx(survival_data_lin, "results/spreadsheets/survival_data.xlsx")
dev.off()
