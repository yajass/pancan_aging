library(wateRmelon)
library(GeneOverlap)
library(pbmcapply)
library(data.table)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(limma)
library(edgeR)
library(TCGAutils)
library(tidyverse)
library(limma)
library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)

# load data
source("code/useful_functions.R")
source("code/02 - survival.R")
survival_data <- openxlsx::read.xlsx("results/spreadsheets/survival_data.xlsx")
conf.data <- openxlsx::read.xlsx("results/spreadsheets/confounder_identification/confounders_to_use.xlsx")
conf.data[is.na(conf.data)] = ""
load("../Data/NewTCGAcounts/Full Matrix.RData")
folder <- "../Data/methylation_gene_beta/"
meth_codes <- gsub("\\..*","", list.files(folder))
types <- unique(survival_data$type)
all(meth_codes == types)

# de functions
deg.linear <- function(cancers = types, counts = RNAseq_counts, confounders = conf.data){
  
  cancer_age_cont <- survival_data %>% dplyr::filter(type == cancers)
  
  # match dataframes
  matched_counts <- matchCountsClinical(filtered_counts = RNAseq_counts, survival_data = cancer_age_cont)
  cancer_counts_cont <- matched_counts[[1]]
  cancer_age_cont <- matched_counts[[2]]
  samples <- matched_counts[[2]]
  rm(counts)
  
  # prep confounders for model matrix
  conf.vars <- confounders$linear[confounders$type == cancers]
  if (nchar(conf.vars) > 2){
    if(stringr::str_detect(conf.vars, pattern = "\\+")){
      conf.vars <- unlist(stringr::str_split(conf.vars, "\\+"))
    }
  }
  if (nchar(conf.vars[1]) < 2){
    conf.vars <- NULL
  }
  cancer_age_cont <- cancer_age_cont %>% dplyr::ungroup() %>%
    dplyr::select(age_at_initial_pathologic_diagnosis, conf.vars)
  
  # run limma on gene expresdsion
  form_lin <- as.formula(paste0("~ ", paste(colnames(cancer_age_cont), collapse = "+"))); print(form_lin)
  design <- model.matrix(form_lin, data = cancer_age_cont)
  cancer_counts_cont <- cancer_counts_cont[which(rowSums(cpm(cancer_counts_cont)> 2) >= 2) ,]
  print(dim(cancer_counts_cont))
  Data.voom <- voom(cancer_counts_cont,plot=FALSE)
  fit <- lmFit(Data.voom,design)
  fit2 <-eBayes(fit)
  deg.res <- topTable(fit2,n=Inf, coef = "age_at_initial_pathologic_diagnosis") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(type = cancers,
                  assay = "rna_seq")
  
  ## methylation
  # load
  meth.dat <- data.table::fread(paste0(folder,list.files(folder, pattern = cancers)))[-1, ] %>%
    tibble::column_to_rownames(., "Hybridization REF")
  dn <- dimnames(meth.dat)
  meth.dat <- apply(meth.dat, 2, as.numeric)
  meth.dat <- data.frame(Beta2M(B = meth.dat))
  
  dimnames(meth.dat) <- dn
  colnames(meth.dat) <- TCGAbarcode(colnames(meth.dat))
  
  # preproc
  common <- intersect(colnames(meth.dat), samples$bcr_patient_barcode)
  meth.dat <- meth.dat %>% dplyr::select(sort(common))
  samples <- samples %>% dplyr::filter(bcr_patient_barcode %in% common) %>%
    dplyr::arrange(bcr_patient_barcode)
  samples_filt <- samples %>% dplyr::ungroup() %>%
    dplyr::select(age_at_initial_pathologic_diagnosis, conf.vars)
  print(all(colnames(meth.dat) == samples$bcr_patient_barcode))
  
  # design mattrix + limma
  print(dim(meth.dat))
  if (cancers != "GBM"){
    design <- model.matrix(form_lin, data = samples_filt)
  }
  if (cancers == "GBM"){
    design <- model.matrix(~age_at_initial_pathologic_diagnosis, data = samples_filt)
  }
  fit <- lmFit(meth.dat, design)
  fit2 <-eBayes(fit)
  dmg.res <- topTable(fit2,n=Inf, coef = "age_at_initial_pathologic_diagnosis") %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(type = cancers,
                  assay = "dna_meth")
  
  return(list(deg = deg.res, dmg = dmg.res))
}

# colllect results
diff.results <- pbmclapply(types, deg.linear, counts=RNASeq_counts, confounders=conf.data, mc.cores = 12)

# save results
diff.res.combined <- Reduce("cbind", diff.results)
dmg.res <- na.omit(setNames(diff.res.combined[seq(0, 64, 2)], types))
deg.res <- na.omit(setNames(diff.res.combined[seq(1,by=2, len=64)], types))[1:32]; rm(diff.res.combined); rm(diff.results)
saveRDS(deg.res, "results/rds/linear/deg_ensembl.rds")
saveRDS(dmg.res, "results/rds/linear/dmg_gene.rds")
deg.res <- lapply(deg.res, function(x) getGeneNamesFromENSELBL(tibble::column_to_rownames(x, "gene")))
saveRDS(deg.res, "results/rds/linear/deg_gene.rds")
