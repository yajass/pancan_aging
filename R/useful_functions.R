# color scheme
cols <- c(setNames(RColorBrewer::brewer.pal('Dark2', n=3)[1:2], c('Not Age-Associated', 'Age-Associated')),
          setNames(RColorBrewer::brewer.pal('Dark2', n=3)[1:2], c('NAA', 'AA')),
          setNames(ggsci::pal_aaas()(2), c('Tumor', 'Normal')),
          setNames(ggsci::pal_jco()(3), c('Young', 'Old', 'NS')),
          setNames(ggsci::pal_jama()(2), c('FDR < 0.05', 'NS')))


# define age associated and not age associated
aa_cancers <- c('BRCA','LGG','LUSC','OV','THCA','UCEC')
naa_cancers <- c('BLCA','COAD','GBM','HNSC','KIRC','LAML','PAAD','KICH','SARC', 'STAD') ## removed STAD
total_cancers <- c(aa_cancers, naa_cancers)

# immunotherapy genes
icb.genes = c('BTLA','CD27','CD274','CD276','CD40','CD40LG','CD70','CTLA4','ENTPD1','FGL1','HAVCR2','HHLA2','ICOS','ICOSLG','IDO1','LAG3','NCR3','NT5E','PDCD1','PDCD1LG2','SIGLEC15','TM1GD2',
              'TNFRSF18','TNFRSF4','TNFRSF9','TNFSF14','VTCN1')

matchCountsClinical <- function(filtered_counts = NULL, survival_data = NULL, onlyPrimary = TRUE, primaryAndNormal = FALSE){
  library(TCGAbiolinks)
  library(TCGAutils)
  
  if (onlyPrimary){
    # only look at primary solid tumors
    filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition %in% c("Primary Solid Tumor",
                                                                                                          "Primary Blood Derived Cancer - Peripheral Blood" ))]
  }
  if (primaryAndNormal){
    # only look at primary solid tumors
    filtered_counts <- filtered_counts[ ,(TCGAbiospec(colnames(filtered_counts))$sample_definition %in% c("Primary Solid Tumor",
                                                                                                          "Primary Blood Derived Cancer - Peripheral Blood",
                                                                                                          "Blood Derived Normal",
                                                                                                          "Solid Tissue Normal",
                                                                                                          "Bone Marrow Normal",
                                                                                                          "Buccal Cell Normal"))]
  }
  # Remove duplicates if they exist
  biospec_table <- TCGAbiospec(colnames(filtered_counts))
  biospec_table <- biospec_table[order(biospec_table$plate,decreasing = TRUE), ]
  filtered_counts <- filtered_counts[,colnames(filtered_counts)[order(biospec_table$plate, decreasing = TRUE)]]
  if (sum(duplicated(biospec_table$submitter_id)) > 0){
    filtered_counts <- filtered_counts[, !duplicated(biospec_table$submitter_id)]
  }
  
  # change column names to match metadata
  colnames(filtered_counts) <- TCGAutils::TCGAbarcode(colnames(filtered_counts))
  filtered_counts <- filtered_counts[, nchar(colnames(filtered_counts)) == 12]
  
  ## selecting case IDs common to both datasets and discarding the rest
  common_cases <- intersect(colnames(filtered_counts), survival_data$bcr_patient_barcode)
  filtered_age <- survival_data[survival_data$bcr_patient_barcode %in% common_cases,]
  filtered_counts <- filtered_counts[,colnames(filtered_counts) %in% common_cases]
  
  ## sorting case IDs
  filtered_age <- filtered_age[order(filtered_age$bcr_patient_barcode),]
  filtered_counts <- filtered_counts[,order(colnames(filtered_counts))]
  
  if(!all(filtered_age$bcr_patient_barcode == colnames(filtered_counts))){
    print("SAMPLES DO NOT MATCH UP!")
  }
  
  return(list(filtered_counts = filtered_counts, survival_data = filtered_age))
}


getGeneNamesFromENSELBL <- function(countData){
  library(org.Hs.eg.db)
  # rownames(countData) <- gsub("\\..*","",rownames(countData))
  countData[,ncol(countData)+1] <- mapIds(org.Hs.eg.db, keys = gsub("\\..*","",rownames(countData)),
                                          keytype = "ENSEMBL", column="SYMBOL", "first")
  countData <- na.omit(countData)
  countData <- countData[-which(duplicated(countData[,ncol(countData)])),]
  rownames(countData) <- countData[,ncol(countData)]
  countData <- countData[,-ncol(countData)]
  return(countData)
}

runGSEA <- function(diff.results, pathways = gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GSEA_pathways/h.all.v7.0.symbols.gmt")){
  library(fgsea)
  library(magrittr)
  # pathways.hallmark <- gmtPathways("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/GSEA_pathways/h.all.v7.0.symbols.gmt")
  set.seed(2019)
  return(fgsea(pathways = pathways,
               stats = setNames(diff.results %>% dplyr::arrange(t) %>% dplyr::pull(t), 
                                diff.results %>% dplyr::arrange(t) %>% dplyr::pull(gene)),
               nperm = 1000))
}

# run a linear model for cell ~ age + confounders, run pearson correlation
cibersortLm <- function(cell_type, complete_df, conf){
  library(tidymodels)
  # setup confounder formula
  conf = conf[conf$type == unique(complete_df$type),]$linear
  if(conf == ""){
    conf = NULL
    form1 <- as.formula("response ~ age_at_initial_pathologic_diagnosis")
  } else {
    form1 <- as.formula(paste0("response ~ age_at_initial_pathologic_diagnosis +", conf))
  }
  
  # get resposne
  response <- complete_df %>% dplyr::select(cell_type) %>% pull()
  
  # make response normally distributed
  response <- RNOmni::rankNorm(response)
  
  # make lm
  # mod <- lm(form1, data = complete_df); plot(mod)
  mod <- linear_reg() %>%
    set_mode("regression") %>%
    set_engine("lm") %>%
    fit(form1, complete_df)
  
  # do correlation
  corr <- cor.test(response, complete_df$age_at_initial_pathologic_diagnosis)
  cor.df <- data.frame(method = "Pearson", cor = corr$estimate, ci.low =corr$conf.int[1],
                       ci.up = corr$conf.int[2], p = corr$p.value, cancer = unique(complete_df$type),
                       cell_type = cell_type,
                       row.names = NULL)
  return(list(lm_res = tidy(mod) %>% 
                dplyr::mutate(cancer = unique(complete_df$type),
                              cell_type = cell_type),
              cor_res = cor.df))
}
