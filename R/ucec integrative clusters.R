library(TCGAbiolinks)
library(ggpubr)
library(maftools)
library(dplyr)

source("code/useful_functions.R")
source("code/02 - survival.R")

# Data from https://www.nature.com/articles/nature12113
# Integrated genomic characterization of endometrial carcinoma
# Integrative cluster based on CNV, RNASeq, mutations and DNA methylation
data <- readxl::read_xls("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA UCEC/datafile.S1.1.KeyClinicalData.xls")

# Get clinical data and quartile
tumor_clinical <- survival_data %>%
  filter(type == "UCEC", Quartiles %in% c(1,4)) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))

# Integrate datasets
young_samples <- tumor_clinical$bcr_patient_barcode[tumor_clinical$Quartiles == "Young"]
old_samples <- tumor_clinical$bcr_patient_barcode[tumor_clinical$Quartiles == "Old"]
data <- data %>%
  mutate(Quartiles = ifelse(bcr_patient_barcode %in% young_samples, "Young",
                            ifelse(bcr_patient_barcode %in% old_samples, "Old", NA))) %>%
  drop_na()
data$IntegrativeCluster

# Dataframe to plot
df <- data.frame(table(data$Quartiles, data$IntegrativeCluster))
df$Var1 <- factor(df$Var1, levels = c("Young", "Old"))

# visualize
ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_fill_brewer() +
  # theme_pubclean(base_size = 35) +
  theme_transparent(40) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right", legend.box = 'vertical') +
  labs(fill = NULL) +
  guides(fill = guide_legend(nrow = 5))
ggsave("results/figures/fig4/uecc_integrative_clusters.eps",
       dpi = 320, height = 15, width = 8)

df %>% reshape2::dcast(Var2~Var1, value.var = "Freq") %>%
  column_to_rownames("Var2") %>%
  rstatix::row_wise_fisher_test(p.adjust.method = "fdr", detailed = TRUE)
