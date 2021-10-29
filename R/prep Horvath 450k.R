library(pbmcapply)
library(data.table); setDTthreads(1)
library(TCGAbiolinks)
library(TCGAutils)

wd = getwd()
source("code/useful_functions.R")
source("code/02 - survival.R")

# get cancers
total_cancers <- paste0("TCGA-",c(aa_cancers, naa_cancers))


#### RUN THIS THE FIRST TIME ONLY!
# download 450k (non-ffpe) data for total_cncers from gdac firehose
# unpack .tar.gz directory ("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/data/")
# will need to rename coadread to coad manually
# parse downloaded folders and prepare - just leave txt/csv files in the data dir
dir <- ("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/data/")
moveFile <- function(directory){
  system(paste('cd',directory,';
               rm MANIFEST.txt;
               mv * ..;
               cd ..;
               rm', directory))
}
renameFile <- function(files){
  system(paste0('mv ', files,' ', str_extract(files, "[^.]+"),'_methylattion_450.txt'))
}
pbmclapply(list.dirs(dir)[2:length(list.dirs(dir))], moveFile, mc.cores = 5)
pbmclapply(list.files(dir), renameFile, mc.cores = 5)


## Adjust the datasets to fit horvath's input
## trim to keep only the ~350 probes
## can do advanced analysis with laml - use whole 450k for t that
datMiniAnnotation = read.csv("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/datMiniAnnotation3.csv")

trimAndSave <- function(files){
  
  code <- stringr::str_extract(files, "[^.]+")
  # read and remove non beta columns
  meth = fread(paste0(dir, files))
  meth <- data.frame(meth)[, c(1, seq(2, ncol(meth)-1, 4))] %>%
    dplyr::rename(ProbeID = Hybridization.REF)
  colnames(meth) <- gsub("\\.","-",colnames(meth))
  
  # trim to match probes
  match1=match(datMiniAnnotation[,1],meth[,1])
  meth[match1,] %>% dplyr::mutate(ProbeID = as.character(ProbeID))
  
  ##### match with survival data - all samples
  meta <- data.frame(SampleID = colnames(meth)[-1], stringsAsFactors = F) %>%
    dplyr::mutate(bcr_patient_barcode = TCGAbiospec(SampleID)$submitter_id)
  age_sub <- survival_data %>% 
    dplyr::select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis,
                  gender, race, clinical_stage, histological_type, histological_grade) %>%
    dplyr::rename(Age = age_at_initial_pathologic_diagnosis)
  
  # write data
  meth = meth[match1,] %>% dplyr::mutate(ProbeID = as.character(ProbeID))
  meta <- merge(meta, age_sub)[,c(2,1,3,4,5,6,7,8)]
  all(colnames(meth)[-1] == meta$SampleID)
  loc <- "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/input/"
  write.csv(meth, file = paste0(loc,code,"_horvath_input_beta.csv"),
            row.names = FALSE, quote = FALSE)
  write.csv(meta, file = paste0(loc,code,"_horvath_input_meta.csv"),
            row.names = FALSE, quote = FALSE)
  
}
dir <- ("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/data/")
pbmclapply(list.files(dir), trimAndSave, mc.cores = 2, mc.preschedule = TRUE)

# above codde did not work properly - modify output to make sure eveything fits horvath webset
dir = "/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/input/"
betas = list.files(dir)[endsWith(list.files(dir), "beta.csv")]
metas = list.files(dir)[endsWith(list.files(dir), "meta.csv")]
dat21 <- fread("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/probeAnnotation21kdatMethUsed.csv")
correctValues <- function(cancer_code){
  beta = data.frame(fread(paste0("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/input/",
                                 cancer_code,"_horvath_input_beta.csv"), fill = TRUE), stringsAsFactors = F)
  meta = data.frame(fread(paste0("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/input/",
                                 cancer_code,"_horvath_input_meta.csv"), fill = TRUE), stringsAsFactors = F)
  colnames(beta)[-1] <- gsub("\\-",".",colnames(beta)[-1])
  meta$SampleID <- gsub("\\-",".",meta$SampleID)
  beta <- beta %>% dplyr::filter(ProbeID %in% dat21$Name)
  common = intersect(colnames(beta)[-1], meta$SampleID)
  # subset common
  beta <- beta[, c(1, which(colnames(beta) %in% common))]
  beta <- beta[, c(1, order(colnames(beta)[-1]) + 1)]
  meta <- meta %>% dplyr::filter(SampleID %in% common) %>% dplyr::arrange(SampleID)
  message(all(colnames(beta)[-1] == meta$SampleID))
  
  write.csv(beta, paste0("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/input/",
                         cancer_code,"_horvath_input_beta.csv"), row.names = FALSE, quote = FALSE)
  write.csv(meta, paste0("/Users/Yajas/Documents/Elemento/tcga_aging_final/Data/methylation_450/horvath_clock/input/",
                         cancer_code,"_horvath_input_meta.csv"), row.names = FALSE, quote = FALSE)
}
pbmclapply(stringr::str_extract(betas, "[^_]+"), correctValues, mc.cores = 10)

