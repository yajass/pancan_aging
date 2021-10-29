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
aa_res$results %>% dplyr::group_by(type) %>%
  dplyr::summarize(total = length(adjPval), 
                   SigYoung = sum(adjPval < 0.05 & or < 1),
                   SigOld = sum(adjPval < 0.05 & or > 1))

## preproc results and plot
res <- aa_res$results %>%
  mutate(combined = paste(Hugo_Symbol, type, sep="-")) %>%
  inner_join(drivers) %>%
  filter(adjPval < 0.05) %>%
  mutate(or = log10(or), ci.up = log10(ci.up), ci.low = log10(ci.low), name = paste(Hugo_Symbol, type, sep='-'),
         dir = factor(ifelse(or < 0, "Young", "Old"), levels = c("Young","Old"))) %>%
  arrange(or)
ucec.genes <- res %>% filter(or < 0, Cancer == 'UCEC') %>% pull(Gene)

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
ggsave("results/figures/fig4/pancan_aa_forest.eps", height = 13, width = 19)

## Pathway analysis of UCEC mutations
ucec_old <- aa_res$results %>%
  mutate(combined = paste(Hugo_Symbol, type, sep="-")) %>%
  inner_join(drivers) %>%
  filter(adjPval < 0.05, or > 1, type == "UCEC") %>%
  pull(Hugo_Symbol)
ucec_young <- aa_res$results %>%
  mutate(combined = paste(Hugo_Symbol, type, sep="-")) %>%
  inner_join(drivers) %>%
  filter(adjPval < 0.05, or < 1, type == "UCEC") %>%
  pull(Hugo_Symbol)

mart <- bitr(ucec_young, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
ucec_young <- mart$ENTREZ
mart <- bitr(ucec_old, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
ucec_old <- mart$ENTREZ
res_young <- enrichPathway(gene = ucec_young, organism = "human")
res_old <- enrichPathway(gene = ucec_old, organism = "human")

df = data.frame(Age = c(rep('Young', length(ucec_young)), rep('Old', length(ucec_old))),
                Gene = c(ucec_young, ucec_old),
                stringsAsFactors = FALSE)
res = clusterProfiler::compareCluster(Gene ~ Age, data = df, fun = 'enrichPathway')
dotplot(res) + scale_color_viridis_c()

enrichplot::pairwise_termsim(res_young) %>% enrichplot::emapplot() + 
  scale_color_viridis_c() + labs(color = 'FDR'); ggsave("results/figures/fig4/young_ucec_pathways.jpeg", width = 7, height = 7)
enrichplot::pairwise_termsim(res_old) %>% enrichplot::emapplot() + 
  scale_color_viridis_c() + labs(color = 'FDR'); ggsave("results/figures/fig4/old_ucec_pathways.jpeg", width = 7, height = 7)


## Pathway analysis of ALL mutations
ucec_old <- aa_res$results %>%
  mutate(combined = paste(Hugo_Symbol, type, sep="-")) %>%
  inner_join(drivers) %>%
  filter(adjPval < 0.05, or > 1) %>%
  pull(Hugo_Symbol)
ucec_young <- aa_res$results %>%
  mutate(combined = paste(Hugo_Symbol, type, sep="-")) %>%
  inner_join(drivers) %>%
  filter(adjPval < 0.05, or < 1) %>%
  pull(Hugo_Symbol)

mart <- bitr(ucec_young, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
ucec_young <- mart$ENTREZ
mart <- bitr(ucec_old, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Hugo_Symbol", "ENTREZ")
ucec_old <- mart$ENTREZ
df = data.frame(Age = c(rep('Young', length(ucec_young)), rep('Old', length(ucec_old))),
                Gene = c(ucec_young, ucec_old),
                stringsAsFactors = FALSE)
res = clusterProfiler::compareCluster(Gene ~ Age, data = df, fun = 'enrichPathway')
y <- as.data.frame(res@compareClusterResult) %>%
  mutate(GeneRatio = as.numeric((str_split(GeneRatio, '/') %>% map(1)))/
           as.numeric((str_split(GeneRatio, '/') %>% map(2))))
y = data.frame(y[c(1:4,10,120,121, 126, 127, 123, 124, 136, 149, 151),]) %>% arrange(Cluster, -GeneRatio, p.adjust) %>%
  mutate(Description = factor(Description, levels = unique(Description)))
ggplot(y, aes(x = Cluster, size = GeneRatio, y = Description, color = p.adjust)) + 
  geom_point() +
  theme_bw(base_size = 24) +
  scale_colour_viridis_c() +
  labs(y = NULL, x = NULL, color = 'FDR') +
  scale_size(limits = c(0.01,0.4), range = c(2, 8)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20))
ggsave("results/figures/fig4/all_snv_clusterprofiler_aa.eps", height = 13, width = 12)


## Following on R 4.0
## UCEC Co-oncoplot and lollipop
load('/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/tcga_mutations/UCEC.RData')
ucec_clin = survival_data %>% filter(type == 'UCEC') %>% dplyr::select(bcr_patient_barcode, Quartiles) %>% dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode)
tcga_ucec_mc3@clinical.data = merge(tcga_ucec_mc3@clinical.data, ucec_clin) %>%
  mutate(Quartiles = ifelse(Quartiles == 1, 'Young', ifelse(Quartiles == 2, 'Old', 'Middle Aged')))
young = subsetMaf(tcga_ucec_mc3, clinQuery = "Quartiles == 'Young'")
old = subsetMaf(tcga_ucec_mc3, clinQuery = "Quartiles == 'Old'")


## high impact UCEC mutations
y.high = subsetMaf(young, genes = ucec.genes, query = 'IMPACT == "HIGH"')
o.high = subsetMaf(old, genes = ucec.genes, query = 'IMPACT == "HIGH"')
comp = maftools::mafCompare(m1 = y.high, m2 = o.high); comp

png("results/figures/fig4/ucec_high_impact_cobarplot.png", res = 100)
maftools::coBarplot(m1 = y.high, m2 = o.high, m1Name = 'Young', m2Name = 'Old',
                    pctSize = 1.5, titleSize = 1, geneSize = 1, geneMar = 4,
                    axisSize = 1.2, legendTxtSize = 1.2,
                    genes = c('CTCF','ATRX','BRCA1','BRCA2','FAT1'), orderBy = 'm1')
dev.off()


png("results/figures/fig4/ucec_high_impact_colollipop.png", res = 200, width = 8, height = 4,units = 'in')
lollipopPlot2(m1 = y.high, m2 = o.high, gene = 'ATRX', roundedRect = T,
              showDomainLabel = F, m1_name = 'Young', m2_name = 'Old', alpha = 0.67)
dev.off()
## end of high impact ucec mutations



res %>% filter(Cancer == 'UCEC') %>% dplyr::arrange(adjPval) %>% pull(Hugo_Symbol) %>% head(8)

coOncoplot(m1 = young, m2 = old, m1Name = 'Young', m2Name = 'Old', genes = ucec.genes)
# lollipopPlot2(m1 = young, m2 = old, gene = "TSC1", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "Young", m2_name = "Old",
#               pointSize = 2, legendTxtSize = 2, showDomainLabel = F, roundedRect = TRUE, m1_label = c('891', '203'))

jpeg('results/figures/fig4/lollipop_ucec_p75.jpeg', width = 6, height = 4.5, units = 'in', res = 320)
lollipopPlot2(m1 = young, m2 = old, gene = "PSIP1", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "Young", m2_name = "Old",
              pointSize = 2.5, legendTxtSize = 2, showDomainLabel = F, roundedRect = TRUE, m1_label = c(344, 132, 314), alpha = 3/4)
dev.off()

postscript('results/figures/fig4/cobarplot_ucec.eps', width = 6, height  = 3)
coBarplot(m1 = young, m2 = old, m1Name = "Young", m2Name = "Old", 
          pctSize = 0.8,axisSize = 1.5, titleSize = 1.5,
          normalize = T, orderBy = 'm1', legendTxtSize =0,
          genes = res %>% filter(Cancer == 'UCEC') %>% dplyr::arrange(adjPval) %>% pull(Hugo_Symbol) %>% head(8))
dev.off()


## Drugs - UCEC young
library(ggalluvial)
cairo_ps('results/figures/fig4/drugs_young_ucec.eps', width = 7, height = 8)
young.drug = drugInteractions(maf = young, fontSize = 0.75, genes = ucec.genes)
dev.off()
young.drug.list = drugInteractions(maf = young, fontSize = 0.75, genes = ucec.genes, drugs = T)

## clinically actionable
dat = young.drug.list %>% filter(PMIDs != '', Gene %in% (young.drug$Gene[young.drug$category == 'CLINICALLY ACTIONABLE'])) %>% 
  dplyr::count(Gene, interaction_claim_source, drug_name) %>%
  filter(drug_name != '')
levs = dat %>% dplyr::count(Gene) %>% arrange(-n) %>% pull(Gene)
dat$Gene = factor(dat$Gene, levels = rev(levs))
dat %>% ggplot(aes(x = Gene, y = n, color = interaction_claim_source, fill = interaction_claim_source)) +
  geom_bar(stat = 'identity') + coord_flip() + labs(x = NULL, y = 'No. of Drugs', fill = 'Source', color = 'Source') +
  theme_pubr(20) + theme(legend.position = c(0.8,0.3)) +
  guides(fill = guide_legend(override.aes = list(size=20)))
ggsave('results/figures/fig4/drugs_young_ucec_by_source_new.eps', width = 12, height = 13)

dat %>% filter(drug_name != '') %>% dplyr::count(drug_name) %>% arrange(-n) %>% slice(1:10) %>%
  ggplot(aes(x = n, y = reorder(drug_name, n))) + geom_col() + theme_pubr(20) +
  labs(x = 'No. of Genes', y = NULL)
ggsave('results/figures/fig4/drugs_young_ucec_by_drug_new.eps', width = 12, height = 13)


## drug categories
dat = young.drug.list %>% filter(PMIDs != '') %>% 
  dplyr::count(Gene, interaction_claim_source, drug_name)
genes = dat %>% dplyr::count(Gene) %>% arrange(-n) %>% pull(Gene)
young.drug %>% filter(Gene %in% genes) %>% dplyr::count(category) %>% arrange(-n) %>%
  ggplot(aes(x = n, y = reorder(category, n))) +
  geom_bar(stat = 'identity') + theme_pubr(20) + 
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  labs(y = NULL, x = 'No. of Genes') + theme_pubr(20)
ggsave('results/figures/fig4/drugs_young_ucec_by_category_new.eps', width = 12, height = 13)


