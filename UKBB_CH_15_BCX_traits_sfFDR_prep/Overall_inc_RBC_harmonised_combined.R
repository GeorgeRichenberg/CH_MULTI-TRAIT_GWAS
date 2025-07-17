### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_rbc <- vroom("overall_inc_RBC_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_hgb <- vroom("overall_inc_HGB_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_hct <- vroom("overall_inc_HCT_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_mch <- vroom("overall_inc_MCH_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_mcv <- vroom("overall_inc_MCV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_mchc <- vroom("overall_inc_MCHC_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_rdw <- vroom("overall_inc_RDW_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

overall_rbc <- rename_bcx_columns(overall_rbc, "rbc")
overall_hgb <- rename_bcx_columns(overall_hgb, "hgb")
overall_hct <- rename_bcx_columns(overall_hct, "hct")
overall_mch <- rename_bcx_columns(overall_mch, "mch")
overall_mcv <- rename_bcx_columns(overall_mcv, "mcv")
overall_mchc <- rename_bcx_columns(overall_mchc, "mchc")
overall_rdw <- rename_bcx_columns(overall_rdw, "rdw")

### Merge all datasets:
overall_bcx <- overall_rbc %>%
  inner_join(overall_hgb, by = "SNP", suffix = c("_rbc", "_hgb")) %>%
  inner_join(overall_hct, by = "SNP", suffix = c("", "_hct")) %>%
  inner_join(overall_mch, by = "SNP", suffix = c("", "_mch")) %>%
  inner_join(overall_mcv, by = "SNP", suffix = c("", "_mcv")) %>%
  inner_join(overall_mchc, by = "SNP", suffix = c("", "_mchc")) %>%
  inner_join(overall_rdw, by = "SNP", suffix = c("", "_rdw"))

### Remove duplicate rows:
overall_bcx <- overall_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- overall_bcx$SNP[duplicated(overall_bcx$SNP)]
duplicates

### Save output:
vroom_write(overall_bcx, file = "overall_inc_bcx_rbc_harmonised_combined.tsv.gz")