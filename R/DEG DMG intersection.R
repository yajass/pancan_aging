library(tidyverse)
library(data.table)
library(GeneOverlap)
library(corrplot)
library(ggcorrplot)
library(patchwork)
library(ggpubr); theme_set(theme_pubr(20))
library(ggsci)
library(ReactomePA)
library(clusterProfiler)

source("code/useful_functions.R")
source("code/02 - survival.R")

## load data ##
dmg <- readRDS("results/rds/logistic/dmg_gene.rds")[aa_cancers]
deg <- lapply(readRDS("results/rds/logistic/deg_gene.rds")[aa_cancers], rownames_to_column, "gene")
cgc = fread('../Data/COSMIC/CGC.tsv', check.names = TRUE)

## Find overlap between tumor types - DEG ##
############################################
# up and down regulated genes
up_genes_rna <- lapply(deg, function(u)u$gene[u$adj.P.Val < 0.05 & u$logFC > 0])
down_genes_rna <- lapply(deg, function(u)u$gene[u$adj.P.Val < 0.05 & u$logFC < 0])

# genes up or down regulated in `i` or more tumor types
i <- 3
old = plyr::ldply(up_genes_rna, cbind,.id = "type") %>% 
  dplyr::rename(gene = `1`) %>%
  dplyr::count(gene) %>%
  dplyr::filter(n >= i) %>% dplyr::arrange(-n) %>%
  mutate(Age = 'Old')
young = plyr::ldply(down_genes_rna, cbind,.id = "type") %>% 
  dplyr::rename(gene = `1`) %>%
  dplyr::count(gene) %>%
  dplyr::filter(n >= i) %>% dplyr::arrange(-n) %>%
  mutate(Age = 'Young')

## check which tumor types the genes are mutated in
library(pbmcapply)
young.gene.types = mclapply(as.character(young$gene), function(this.gene){
  this.gene = as.character(this.gene)
  BRCA = sum(this.gene %in% down_genes_rna$BRCA)
  LGG = sum(this.gene %in% down_genes_rna$LGG)
  LUSC = sum(this.gene %in% down_genes_rna$LUSC)
  OV = sum(this.gene %in% down_genes_rna$OV)
  THCA = sum(this.gene %in% down_genes_rna$THCA)
  UCEC = sum(this.gene %in% down_genes_rna$UCEC)
  return(data.table(gene = as.character(this.gene), BRCA = BRCA, LGG = LGG, LUSC = LUSC, OV = OV, THCA = THCA, UCEC = UCEC))
}, mc.cores = 8) %>% rbindlist
old.gene.types = mclapply(as.character(old$gene), function(this.gene){
  this.gene = as.character(this.gene)
  BRCA = sum(this.gene %in% up_genes_rna$BRCA)
  LGG = sum(this.gene %in% up_genes_rna$LGG)
  LUSC = sum(this.gene %in% up_genes_rna$LUSC)
  OV = sum(this.gene %in% up_genes_rna$OV)
  THCA = sum(this.gene %in% up_genes_rna$THCA)
  UCEC = sum(this.gene %in% up_genes_rna$UCEC)
  return(data.table(gene = as.character(this.gene), BRCA = BRCA, LGG = LGG, LUSC = LUSC, OV = OV, THCA = THCA, UCEC = UCEC))
}, mc.cores = 8) %>% rbindlist

young = merge(young, young.gene.types); dim(young)
old = merge(old, old.gene.types); dim(old)

cgc.age = rbind(merge(young, cgc, by.x = 'gene', by.y = 'Gene.Symbol'),
                merge(old, cgc, by.x = 'gene', by.y = 'Gene.Symbol'))
cgc.age %>% dplyr::select(gene, Age, BRCA, LGG, LUSC, OV, THCA, UCEC, Tier, Hallmark, Molecular.Genetics, Role.in.Cancer) %>%
  fwrite('results/tables/cgc_genes.tsv', row.names = F, quote = F, sep = '\t')
openxlsx::write.xlsx(cgc.age, 'results/spreadsheets/common_deg_cgc.xlsx', asTable = T)

i <- 4
old = plyr::ldply(up_genes_rna, cbind,.id = "type") %>% 
  dplyr::rename(gene = `1`) %>%
  dplyr::count(gene) %>%
  dplyr::filter(n >= i) %>% dplyr::arrange(-n) %>%
  mutate(Age = 'Old')
young = plyr::ldply(down_genes_rna, cbind,.id = "type") %>% 
  dplyr::rename(gene = `1`) %>%
  dplyr::count(gene) %>%
  dplyr::filter(n >= i) %>% dplyr::arrange(-n) %>%
  mutate(Age = 'Young')


# function to create gene overlap plot - pairwise
createOverlap <- function(gene_list, padj_method = "fdr"){
  library(GeneOverlap)
  or_data <- log10(getMatrix(newGOM(gene_list), name = "odds.ratio")+1)
  pval_before_adjust <- getMatrix(newGOM(gene_list),name = "p")
  fdr_data <- matrix(p.adjust(pval_before_adjust, padj_method), 
                     nrow = nrow(pval_before_adjust), 
                     ncol = ncol(pval_before_adjust))
  dimnames(fdr_data) <- dimnames(pval_before_adjust)
  return(list(OR = or_data, fdr = fdr_data))
}

# get overlaps
or_rna_up <- createOverlap(up_genes_rna)$OR; fdr_rna_up <- createOverlap(up_genes_rna)$fdr
or_rna_down <- createOverlap(down_genes_rna)$OR; fdr_rna_down <- createOverlap(down_genes_rna)$fdr

# plot overlaps
corrplot(or_rna_up, is.corr = F, method = "color", p.mat = fdr_rna_up,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.col = "white", tl.col = "black",
         type = "upper", tl.cex = 1.5, outline = "white", col = RColorBrewer::brewer.pal(9, "Blues"))
corrplot(or_rna_down, is.corr = F, method = "color", p.mat = fdr_rna_down,
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.col = "white", tl.col = "black",
         type = "upper", tl.cex = 1.5, outline = "white", col = RColorBrewer::brewer.pal(9, "Reds"))

ggcorrplot(or_rna_up, p.mat = fdr_rna_up, type = 'upper', lab = TRUE, outline.color = 'white', pch.cex = 20) + 
  scale_fill_gradient(low = 'ivory', high = 'red') + labs(fill = 'log OR', title = 'Older') + theme(legend.position = 'top') +
ggcorrplot(or_rna_down, p.mat = fdr_rna_down, type = 'upper', lab = TRUE, outline.color = 'white', pch.cex = 20) +
  scale_fill_gradient(low = 'ivory', high = 'blue') + labs(fill = 'log OR', title = 'Younger') + theme(legend.position = 'top') +
  plot_layout(ncol = 1)
ggsave('results/figures/fig2/rna_overlap_or.eps', width = 7, height = 9)
dev.off()

## Gene overlap - RNASeq and DNAm ##
####################################
up_genes_rna <- lapply(deg, function(u)u$gene[u$adj.P.Val < 0.05 & u$logFC > 0])
down_genes_rna <- lapply(deg, function(u)u$gene[u$adj.P.Val < 0.05 & u$logFC < 0])
hyper_meth_genes <- lapply(dmg, function(u)u$gene[u$adj.P.Val < 0.05 & u$logFC > 0])
hypo_meth_genes <- lapply(dmg, function(u)u$gene[u$adj.P.Val < 0.05 & u$logFC < 0])

names(up_genes_rna) <- paste0(names(up_genes_rna), " mRNA")
names(down_genes_rna) <- paste0(names(down_genes_rna), " mRNA")
names(hyper_meth_genes) <- paste0(names(hyper_meth_genes), " DNAm")
names(hypo_meth_genes) <- paste0(names(hypo_meth_genes), " DNAm")

gom.object_up <- newGOM(lapply(up_genes_rna, unlist),lapply(hypo_meth_genes, unlist))
or_data_up <- log10(getMatrix(gom.object_up, name = "odds.ratio")+1)
pval_data_up_before_adjust <- getMatrix(gom.object_up,name = "p")
pval_data_up <- matrix(p.adjust(pval_data_up_before_adjust, "fdr"), nrow = nrow(pval_data_up_before_adjust), ncol = ncol(pval_data_up_before_adjust))
dimnames(pval_data_up) <- dimnames(pval_data_up_before_adjust)

gom.object_down <- newGOM(lapply(down_genes_rna, unlist),lapply(hyper_meth_genes, unlist))
or_data_down <- log10(getMatrix(gom.object_down, name = "odds.ratio")+1)
pval_data_down_before_adjust <- getMatrix(gom.object_down,name = "p")
pval_data_down <- matrix(p.adjust(pval_data_down_before_adjust, "fdr"), nrow = nrow(pval_data_down_before_adjust), ncol = ncol(pval_data_down_before_adjust))
dimnames(pval_data_down) <- dimnames(pval_data_up_before_adjust)

newDF <- data.frame(OR = c(diag(or_data_down), diag(or_data_up)),
                    FDR = c(p.adjust(diag(pval_data_down_before_adjust), "fdr"),
                            p.adjust(diag(pval_data_up_before_adjust), "fdr")),
                    Tumor = rep(gsub('.{5}$', '', colnames(or_data_down)), 2),
                    Type = c(rep("Young", length(aa_cancers)), rep("Old", length(aa_cancers)))) %>%
  dplyr::mutate(Sig = ifelse(FDR < 0.05, "< 0.05", "≥ 0.05"),
                Type = factor(Type, levels = c("Young", "Old")))
ggbarplot(data = newDF, x = 'Tumor', y = 'OR', fill = 'Sig',
          palette = 'jama') +
  facet_wrap(vars(Type), nrow = 2) +
  labs(x = NULL, y = expression(log[10]~(OR+1)), fill = "FDR") +
  theme_pubr() + coord_flip() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/fig2/deg_dmg_overlap_or.jpg", height = 8, width = 3)
dev.off()


# Plot differential methylation
########################################################################################################################
# mutate results
results_meth_genes <- bind_rows(dmg, .id = "CancerCode")
results_meth_genes <- results_meth_genes %>%
  filter(adj.P.Val < 0.05) %>%
  # filter(abs(logFC) > 0.5) %>%
  group_by(CancerCode) %>%
  mutate(Status = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated")) %>%
  mutate(Case = paste0(CancerCode,"-",Status)) %>%
  mutate(Count = 1)
results_rnaseq_genes <- bind_rows(deg, .id = "CancerCode")
results_rnaseq_genes <- results_rnaseq_genes %>%
  filter(adj.P.Val < 0.05) %>%
  # filter(abs(logFC) > 0.5) %>%
  group_by(CancerCode) %>%
  mutate(Status = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
  mutate(Case = paste0(CancerCode,"-",Status)) %>%
  mutate(Count = 1)

df1 <- as.data.frame(table(results_meth_genes$CancerCode))
df1 <- df1[order(df1$Freq), ]
plotOrder <- factor(df1$Var1, levels = df1$Var1)
plot1 <- ggbarplot(data = df1, x = 'Var1', y = 'Freq', order = plotOrder, fill = "black") +
  labs(x = NULL, y = "DMG") +
  theme_pubclean(base_size = 35) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank())

df2 <- as.data.frame(table(results_meth_genes$CancerCode, results_meth_genes$Status))
df2$Var2 <- factor(df2$Var2, levels = c("Hypomethylated","Hypermethylated"))
plot2 <- ggbarplot(data = df2, x = 'Var1', y = 'Freq', fill = 'Var2', position = position_fill(),
                   order = plotOrder) +
  labs(x = NULL, y = "Fraction of DMG", fill = "Methylation") +
  ggsci::scale_fill_jco(labels = c("Young","Old")) + 
  theme_pubclean(base_size = 35) +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank())
ggarrange(plot1,plot2, ncol = 1, align = "h", common.legend = T, heights = c(1.3,2))
ggsave("results/figures/fig2/dmg_fractions.eps")


# clusterProfiler Pathway Analysis
########################################################################################################################
# combine meth and RNASeq results
colnames(results_meth_genes)[3:ncol(results_meth_genes)] <- paste0(colnames(results_meth_genes)[3:ncol(results_meth_genes)], "_DNAm")
combined_results <- results_meth_genes %>% right_join(results_rnaseq_genes, by = c("CancerCode", "gene"))
combined_results <- na.omit(data.frame(CancerCode = combined_results$CancerCode,
                                       Gene = combined_results$gene, 
                                       DNAm = combined_results$Status_DNAm,
                                       Expression = combined_results$Status))
dna.rna.interaction <- combined_results %>%
  dplyr::group_by(CancerCode) %>%
  dplyr::count(DNAm, Expression)

combined_results$Intersection <- paste0(combined_results$DNAm,"-",combined_results$Expression)
accepted <- c("Hypermethylated-Downregulated", "Hypomethylated-Upregulated")
combined_results <- combined_results[combined_results$Intersection %in% accepted, ]
dna.rna.interaction <- combined_results %>%
  dplyr::group_by(CancerCode) %>%
  dplyr::count(DNAm, Expression)

# Genes up/downregulated across cancers
combined_results %>% dplyr::filter(Intersection == "Hypermethylated-Downregulated") %>% dplyr::count(Gene) %>% dplyr::arrange(-n)
combined_results %>% dplyr::filter(Intersection == "Hypomethylated-Upregulated") %>% dplyr::count(Gene) %>% dplyr::arrange(-n)

## Cluster Profiler
combined_results$CancerCode.Intersection <- paste0(combined_results$CancerCode, "-", combined_results$Intersection)
combined_results$CancerCode.Intersection <- as.character(combined_results$CancerCode.Intersection)
keys <- unique(combined_results$CancerCode.Intersection)

# convert to entrez
mart <- bitr(combined_results$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(mart) <- c("Gene", "ENTREZ")
combined_results <- merge(combined_results, mart, "Gene")

gene_list <- setNames(rep(list(NA), length(keys)), keys)
for (i in 1:length(keys)) {
  gene_list[[i]] <- as.character(combined_results$ENTREZ[combined_results$CancerCode.Intersection == keys[i]])
}
names(gene_list) <- gsub('\\-',' ', names(gene_list))

# Overall young vs old clusters
combined_results$Age <- ifelse(combined_results$Expression == "Upregulated", "Old", "Young")
formula_res <- compareCluster(ENTREZ~Age, data=combined_results, fun="enrichPathway")
dotplot(formula_res) + 
  scale_x_discrete(labels = scales::wrap_format(10)) +
  scale_y_discrete(labels = scales::wrap_format(40)) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(5,15)) +
  theme_pubclean(base_size = 25) +
  theme(legend.position = "right", 
        axis.ticks = element_blank())

# Young vs old by tumor type
formula_res <- compareCluster(ENTREZ~CancerCode+Age, data=combined_results, fun="enrichPathway")
strings <- c("senescence", "aging", "immune", "histone", "meth", "telomer", "interferon",
             "clock","il","golgi","repair")
formula_res@compareClusterResult %>%
  fwrite('results/tables/deg_dmg_pathways.tsv', row.names = F, quote = F, sep = '\t')
formula_res@compareClusterResult %>%
  dplyr::filter(grepl(paste(strings,collapse = "|"), Description, ignore.case = TRUE)) %>%
  dplyr::mutate(GeneRatio = round(as.numeric(sub("\\/.*", "", GeneRatio))/as.numeric(sub('.*/', '', GeneRatio)), 2),
                Description = factor(Description, levels = unique(Description))) %>%
  ggplot(., aes(x=Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_size_continuous(range = c(5,15)) +
  scale_x_discrete(labels = scales::wrap_format(4)) +
  scale_y_discrete(labels = scales::wrap_format(40)) +
  ylab(NULL) + xlab(NULL) + scale_color_viridis_c() + labs(color = "FDR") + theme_bw(20)
ggsave("results/figures/fig2/deg_dmg_clusterprofiler_reactome.eps", width = 12, height = 8)
