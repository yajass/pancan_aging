library(clusterProfiler)
library(ReactomePA)
library(pbmcapply)
library(gtable)
library(grid)
library(gtools)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(TCGAbiolinks)
library(maftools)
options(stringsAsFactors = FALSE)

## load data
source("code/useful_functions.R")
source("code/02 - survival.R")
aa_maf <- readRDS("results/rds/try.aa_maf.rds")
naa_maf <- readRDS("results/rds/naa_maf.rds")
drivers <- readxl::read_xlsx("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/TCGA Drivers/driver.genes.xlsx",
                             sheet = 2, trim_ws = TRUE, skip = 3)[,1:2] %>% 
  mutate(Cancer = ifelse(Cancer %in% c("COADREAD", "COAD","READ"), "COAD", Cancer))
drivers <- data.frame(Gene = rep(drivers %>% filter(Cancer == "PANCAN") %>% pull(Gene), each = length(cancer_code)),
                      Cancer = cancer_code) %>%
  bind_rows(drivers) %>%
  filter(Cancer != "PANCAN") %>%
  distinct() %>%
  mutate(combined = paste(Gene,Cancer, sep="-"))

## setup functions 
doStatistics <- function(maf){
  
  cancer_code = as.character(unique(maf@clinical.data$type))
  
  compareMaf <- function(maf, cancer_type){
    
    # subset mafs
    young_maf <- subsetMaf(maf, tsb = maf@clinical.data %>% filter(Quartiles == "Young", type == cancer_type) %>% pull(Tumor_Sample_Barcode) %>% as.character())
    old_maf <- subsetMaf(maf, tsb = maf@clinical.data %>% filter(Quartiles == "Old", type == cancer_type) %>% pull(Tumor_Sample_Barcode) %>% as.character())
    
    # compare young vs old
    possibleError <- tryCatch(
      old_vs_young <- mafCompare(m1 = old_maf, m2 = young_maf, m1Name = 'Old', m2Name = 'Young', minMut = 5),
      error=function(e) e
    )
    if(inherits(possibleError, "error")) return(NULL)
    
    # return results
    return(list(results = data.frame(old_vs_young$results) %>% mutate(type = cancer_type),
                SampleSummary = data.frame(old_vs_young$SampleSummary) %>% mutate(type = cancer_type)))
    
  }
  compared = pbmclapply(cancer_code, compareMaf, maf = maf, mc.cores = 10)
  return(list(results = lapply(compared, function(x) x$results)%>% bind_rows(),
              SampleSummary = lapply(compared, function(x) x$SampleSummary)%>% bind_rows()))
}

## run analysis
aa_res <- doStatistics(aa_maf)
naa_res <- doStatistics(naa_maf)

# overview of results
naa_res$results %>% dplyr::group_by(type) %>%
  dplyr::summarize(total = length(adjPval), 
                   SigYoung = sum(adjPval < 0.05 & or < 1),
                   SigOld = sum(adjPval < 0.05 & or > 1))

## preproc results and plot
res <- naa_res$results %>%
  mutate(combined = paste(Hugo_Symbol, type, sep="-")) %>%
  inner_join(drivers) %>%
  filter(adjPval < 0.05) %>%
  mutate(or = log10(or), ci.up = log10(ci.up), ci.low = log10(ci.low), name = paste(Hugo_Symbol, type, sep='-'),
         dir = factor(ifelse(or < 0, "Young", "Old"), levels = c("Young","Old"))) %>%
  arrange(or)

## plot
forest <- ggplot(res, aes(x = or, xmin = ci.low, xmax = ci.up, y = reorder(name, or), color = dir)) +
  geom_point(size = 6, shape = 18) + geom_errorbar(size = 2) + geom_vline(xintercept = 0, lty = 2) +
  theme_pubclean(30) + scale_color_jco() + xlim(-3,3) +
  labs(x = "log Odds Ratio", y = NULL, color = NULL)
# tables
tab_base <- ggplot(res, aes(y=reorder(name, or))) +
  ylab(NULL) + xlab("  ") +
  theme(plot.title = element_text(hjust = 0.5, size=28), ## centering title on text
        axis.text.x=element_text(color="white"), ## need text to be printed so it stays aligned with figure but white so it's invisible
        axis.line=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
tab1 <- tab_base + 
  geom_text(aes(x=1, label=Young), size = 10) + 
  ggtitle("Young")
tab2 <- tab_base + 
  geom_text(aes(x=1, label=Old), size = 10) + 
  ggtitle("Old")
tab3 <- tab_base + 
  geom_text(aes(x=1, label=gtools::stars.pval(adjPval)), size = 10) + 
  ggtitle("FDR")
ggarrange(forest, tab1, tab2, tab3, widths = c(10,1,1,1), nrow = 1, align = "h",  common.legend = T, legend = F)
ggsave("results/figures/fig4/pancan_naa_forest.eps", height = 6, width = 19)

