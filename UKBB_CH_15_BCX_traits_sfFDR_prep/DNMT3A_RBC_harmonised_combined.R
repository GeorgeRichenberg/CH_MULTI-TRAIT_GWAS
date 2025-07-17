### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
dnmt3a_rbc <- vroom("dnmt3a_RBC_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_hgb <- vroom("dnmt3a_HGB_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_hct <- vroom("dnmt3a_HCT_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_mch <- vroom("dnmt3a_MCH_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_mcv <- vroom("dnmt3a_MCV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_mchc <- vroom("dnmt3a_MCHC_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_rdw <- vroom("dnmt3a_RDW_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

dnmt3a_rbc <- rename_bcx_columns(dnmt3a_rbc, "rbc")
dnmt3a_hgb <- rename_bcx_columns(dnmt3a_hgb, "hgb")
dnmt3a_hct <- rename_bcx_columns(dnmt3a_hct, "hct")
dnmt3a_mch <- rename_bcx_columns(dnmt3a_mch, "mch")
dnmt3a_mcv <- rename_bcx_columns(dnmt3a_mcv, "mcv")
dnmt3a_mchc <- rename_bcx_columns(dnmt3a_mchc, "mchc")
dnmt3a_rdw <- rename_bcx_columns(dnmt3a_rdw, "rdw")

### Merge all datasets:
dnmt3a_bcx <- dnmt3a_rbc %>%
  inner_join(dnmt3a_hgb, by = "SNP", suffix = c("_rbc", "_hgb")) %>%
  inner_join(dnmt3a_hct, by = "SNP", suffix = c("", "_hct")) %>%
  inner_join(dnmt3a_mch, by = "SNP", suffix = c("", "_mch")) %>%
  inner_join(dnmt3a_mcv, by = "SNP", suffix = c("", "_mcv")) %>%
  inner_join(dnmt3a_mchc, by = "SNP", suffix = c("", "_mchc")) %>%
  inner_join(dnmt3a_rdw, by = "SNP", suffix = c("", "_rdw"))

### Remove duplicate rows:
dnmt3a_bcx <- dnmt3a_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- dnmt3a_bcx$SNP[duplicated(dnmt3a_bcx$SNP)]
duplicates

### Save output:
vroom_write(dnmt3a_bcx, file = "dnmt3a_bcx_rbc_harmonised_combined.tsv.gz")