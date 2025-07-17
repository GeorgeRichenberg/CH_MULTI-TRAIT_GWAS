### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
asxl1_rbc <- vroom("asxl1_RBC_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_hgb <- vroom("asxl1_HGB_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_hct <- vroom("asxl1_HCT_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_mch <- vroom("asxl1_MCH_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_mcv <- vroom("asxl1_MCV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_mchc <- vroom("asxl1_MCHC_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_rdw <- vroom("asxl1_RDW_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

asxl1_rbc <- rename_bcx_columns(asxl1_rbc, "rbc")
asxl1_hgb <- rename_bcx_columns(asxl1_hgb, "hgb")
asxl1_hct <- rename_bcx_columns(asxl1_hct, "hct")
asxl1_mch <- rename_bcx_columns(asxl1_mch, "mch")
asxl1_mcv <- rename_bcx_columns(asxl1_mcv, "mcv")
asxl1_mchc <- rename_bcx_columns(asxl1_mchc, "mchc")
asxl1_rdw <- rename_bcx_columns(asxl1_rdw, "rdw")

### Merge all datasets:
asxl1_bcx <- asxl1_rbc %>%
  inner_join(asxl1_hgb, by = "SNP", suffix = c("_rbc", "_hgb")) %>%
  inner_join(asxl1_hct, by = "SNP", suffix = c("", "_hct")) %>%
  inner_join(asxl1_mch, by = "SNP", suffix = c("", "_mch")) %>%
  inner_join(asxl1_mcv, by = "SNP", suffix = c("", "_mcv")) %>%
  inner_join(asxl1_mchc, by = "SNP", suffix = c("", "_mchc")) %>%
  inner_join(asxl1_rdw, by = "SNP", suffix = c("", "_rdw"))

### Remove duplicate rows:
asxl1_bcx <- asxl1_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- asxl1_bcx$SNP[duplicated(asxl1_bcx$SNP)]
duplicates

### Save output:
vroom_write(asxl1_bcx, file = "asxl1_bcx_rbc_harmonised_combined.tsv.gz")