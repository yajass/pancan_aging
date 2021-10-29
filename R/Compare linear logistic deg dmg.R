library(ggpubr); library(ggsci); library(gginnards)
library(data.table); library(tidyverse)

# load data
source("code/useful_functions.R")
source("code/02 - survival.R")
deg_linear <- readRDS("~/Documents/Elemento/tcga_aging_final/logistic_reviwer_response_1/results/rds/linear/deg_ensembl.rds") %>% 
  bind_rows() %>% mutate(model = "linear") %>% filter(type %in% sig_cancers)
deg_logistic <- readRDS("~/Documents/Elemento/tcga_aging_final/logistic_reviwer_response_1/results/rds/logistic/deg_ensembl.rds") %>% 
  bind_rows() %>% mutate(model = "logistic") %>% filter(type %in% sig_cancers)

# deg_linear %>%
#   dplyr::mutate(Sig = ifelse(adj.P.Val < 0.05, "Sig", "NS")) %>%
#   dplyr::group_by(type) %>%
#   dplyr::count(Sig) %>%
#   dplyr::filter(Sig == "Sig") %>%
#   dplyr::bind_rows(data.frame(type = setdiff(sig_cancers, .$type),
#                               Sig = "Sig",
#                               n = 0)) %>%
#   dplyr::inner_join(setNames(data.frame(table(survival_data$type)), c("type", "samp"))) %>%
#   dplyr::mutate(n = 100*n/56537,
#                 age_ass = ifelse(n > 1, "Age-Associated", "Not Age-Associated")) %>%
#   ggplot(., aes(x = reorder(type, n), y = n, label = round(n,2), fill = age_ass)) +
#   geom_bar(stat = "identity") +
#   labs(x = NULL, y = "Percent DEG", fill = "Molecular\nPhenotype") +
#   geom_hline(yintercept = 1 , lty = 2) +
#   scale_fill_aaas() +
#   coord_flip() +
#   theme(legend.position = c(0.6,0.3))

# count sig DEG
deg_linear <- deg_linear %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::count(type) %>%
  dplyr::rename(`Continuous Age` = n)
deg_logistic <- deg_logistic %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::count(type) %>%
  dplyr::rename(`Binned Age` = n)
deg <- merge(deg_logistic, deg_linear, all.x = TRUE) %>%
  melt('type') %>%
  mutate(value = 100*value/56537,
         sel = ifelse(value >= 1, "Age-Associated", "Not Age-Associated"))


p <- ggviolin(deg %>% filter(sel == "Age-Associated") %>% mutate(variable = ifelse(variable == "Binned Age", "Quartiled Age", "Continuous Age")),
              x = "variable", y = "value", fill = "variable", add = "jitter", palette = 'aaas') +
  theme_bw(20) + ylim(-8,25) +  labs(x = NULL, y = "Percent DEG") +
  stat_compare_means(comparisons = list(c('Quartiled Age', 'Continuous Age')), label = "p.signif", label.y = 25, size = 10, bracket.size = 1.5) +
  theme_pubr(20) + theme(legend.position = 'none')
p$layers[[which_layers(p, "GeomSignif")]]$aes_params$textsize <- 10; p
ggsave('results/figures/fig1/binned_vs_continuous_aa_age_percent_deg_violin.eps', width = 8, height = 8)

linear.aa = na.omit(deg$type[deg$sel == "Age-Associated" & deg$variable == "Continuous Age"])
logistic.aa = na.omit(deg$type[deg$sel == "Age-Associated" & deg$variable == "Binned Age"])
linear.naa = na.omit(deg$type[deg$sel == "Not Age-Associated" & deg$variable == "Continuous Age"])
logistic.naa = na.omit(deg$type[deg$sel == "Not Age-Associated" & deg$variable == "Binned Age"])

ggvenn::ggvenn(data = list(`Quartiled Age\nAge-Associated` = logistic.aa, `Continuous Age\nAge-Associated` = linear.aa), 
               fill_alpha = 2/3, stroke_size = 0, show_elements = T, label_sep = "\n", fill_color = pal_aaas()(2))
ggsave('results/figures/fig1/binned_vs_continuous_age_aa_venn.jpeg', width = 6, height = 6)
