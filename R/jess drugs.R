library(tidyverse)
library(RColorBrewer)
library(VennDiagram)
library(ComplexHeatmap)
library(pbmcapply)

#####################
###   Load Data   ###
#####################

data_path <- '/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/jess_results/'

## load MOA dataframe - needed to annotate MOA to output csv
df_new_moa <- read.csv(paste(data_path, 'new_moa.csv', sep = '/'),
                       sep = ",", header = TRUE, stringsAsFactors = FALSE)

## load old and young GSEA sets, total counts for drugs and moas
load(file = paste(data_path, 'yajas_total.Rda', sep = '/'),
     verbose = TRUE)


## filter by p-val and adjusted p-val for down regulated list
create_top_down <- function(input, padj = 0.25, pval = 0.05){
  top_down <- input %>% 
    filter(w < 0, padj_down <= padj, pval_down <= pval) %>% 
    as.data.frame()
  
  top_down <- merge(top_down, df_new_moa[, c("Drugs", "moa_final")], 
                    by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  colnames(top_down)[11] <- "MOA"
  
  #top_down$MOA <- as.character(top_down$MOA)
  top_down$MOA <- unlist(lapply(top_down$MOA, function(x) strsplit(x, ",")[[1]][1]))
  top_down <- top_down[order(top_down$w), ]
  top_down[, c("MOA", "Drug")]
  return(top_down)
}

## filter by p-val and adjusted p-val for both up and down regulated list
create_top_total <- function(input, 
                             down_padj = 0.25, down_pval = 0.05,
                             up_padj = 0.25, up_pval = 0.05){
  top_total <- input %>% 
    filter(w < 0, 
           padj_down <= down_padj, pval_down <= down_pval, 
           padj_up <= up_padj, pval_up <= up_pval) %>% 
    as.data.frame()
  
  top_total <- merge(top_total, df_new_moa[, c("Drugs", "moa_final")], 
                     by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  colnames(top_total)[11] <- "MOA"
  
  #top_total$MOA <- as.character(top_total$MOA)
  top_total$MOA <- unlist(lapply(top_total$MOA, function(x) strsplit(x, ",")[[1]][1]))
  top_total <- top_total[order(top_total$w), ]
  top_total[, c("MOA", "Drug")]
  return(top_total)
}

## filter by p-val and adjusted p-val for down regulated list
## filter by weighted connectivity score 
## don't recommend including up regulated p-vals (will get few, if any, hits)
create_top_weighted <- function(input, 
                                down_padj = 0.25, down_pval = 0.05,
                                # up_padj = 0.25, up_pval = 0.05, 
                                w_max = -0.8){
  top_total <- input %>% 
    filter(w < 0, 
           padj_down <= down_padj, pval_down <= down_pval, 
           # padj_up <= up_padj, pval_up <= up_pval, 
           w <= w_max) %>% 
    as.data.frame()
  
  top_total <- merge(top_total, df_new_moa[, c("Drugs", "moa_final")], 
                     by.x = "Drug", by.y = "Drugs", all.x = TRUE)
  colnames(top_total)[11] <- "MOA"
  
  #top_total$MOA <- as.character(top_total$MOA)
  top_total$MOA <- unlist(lapply(top_total$MOA, function(x) strsplit(x, ",")[[1]][1]))
  top_total <- top_total[order(top_total$w), ]
  top_total[, c("MOA", "Drug")]
  return(top_total)
}

## function to create dataframe of counts
create_counts <- function(df_input, col_name, output_col){
  df_temp <- data.frame(table(df_input[col_name]), stringsAsFactors = FALSE)
  colnames(df_temp) <- c(output_col, 'count')
  df_temp <- df_temp[order(df_temp$count, decreasing = TRUE), ]
  
  ## remove unknown and NA rows before merge (relevant for MOA)
  idx_unknown <- which(df_temp[output_col] == 'Unknown')
  idx_na <- which(is.na(df_temp[output_col]))
  idx_drop <- c(idx_unknown, idx_na)
  if (length(idx_drop)>0){
    return(df_temp[-idx_drop, ])
  } else {
    return(df_temp)
  }
}



####################################################################

list_young_total = lapply(list_fgsea_young, create_top_down)
lapply(list_young_total, nrow)
list_young_drug_count = lapply(list_young_total, function(x) create_counts(x, 'Drug', 'drug'))
lapply(list_young_drug_count, nrow)
list_young_moa_count = lapply(list_young_total, function(x) create_counts(x, 'MOA', 'moa'))
lapply(list_young_moa_count, nrow)

# statistics
fisher_moa_young = pbmclapply(list_young_moa_count, function(x){
  dat = merge(x, df_moa_counts, by = 'moa') %>%
    column_to_rownames('moa')
  colnames(dat) = c('sig_cat', 'total_cat')
  res = rstatix::row_wise_fisher_test(dat, detailed = TRUE, p.adjust.method = 'fdr') %>% arrange(p.adj)
}, mc.cores = 6) %>% bind_rows(.id = 'type')
fisher_drug_young = pbmclapply(list_young_drug_count, function(x){
  dat = merge(x, df_drug_counts, by = 'drug') %>%
    column_to_rownames('drug')
  colnames(dat) = c('sig_cat', 'total_cat')
  res = rstatix::row_wise_fisher_test(dat, detailed = TRUE, p.adjust.method = 'fdr') %>% arrange(p.adj)
}, mc.cores = 6) %>% bind_rows(.id = 'type')


thresh = 0.05
thresh_str = paste0('FDR < ', thresh)
fisher_drug_young %>% filter(p.adj < thresh, estimate > 1) %>% dplyr::count(group) %>% arrange(-n)
drugs = fisher_drug_young %>% filter(p.adj < thresh, estimate > 1) %>% dplyr::count(group) %>% arrange(-n) %>% top_n(5) %>% pull(group)

dat = fisher_drug_young %>% filter(p.adj < thresh, estimate > 1) %>%
  reshape2::dcast(type~group, value.var = 'estimate') %>% column_to_rownames('type')
dat[!is.na(dat)] = thresh_str

# prep for heatmap
ix = order(colSums(dat != thresh_str, na.rm = TRUE))
ix1 = order(rowSums(dat == thresh_str, na.rm = TRUE))
dat = dat[ix1, ix, drop = FALSE]
idx = which(colnames(dat) %in% drugs)

drug_names = HeatmapAnnotation(foo = anno_mark(at = c(idx), labels = colnames(dat)[idx], side = 'bottom'))
anno = df_new_moa %>%
  filter(Drugs %in% colnames(dat)) %>%
  dplyr::select(Drugs, moa_final, Phase) %>%
  column_to_rownames('Drugs') %>%
  dplyr::rename(Mechanism = moa_final) %>%
  mutate(Mechanism = ifelse(str_detect(Mechanism, '\\,'), 'Multiple Targets', Mechanism))
an = HeatmapAnnotation(df = anno, na_col = 'ivory', annotation_name_side = 'left')

pdf('results/figures/fig3/jess_heatmap.pdf', width = 15/1.2, height = 9/1.2)
ComplexHeatmap::oncoPrint(dat, row_labels = toupper(rownames(dat)), bottom_annotation = drug_names, top_annotation = an,
                          width = unit(15, "cm"), height = unit(9, "cm"), show_heatmap_legend = FALSE, use_raster = TRUE, raster_quality = 10)
dev.off()
