### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
dnmt3a_wbc <- vroom("dnmt3a_WBC_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_neu <- vroom("dnmt3a_NEU_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_lym <- vroom("dnmt3a_LYM_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_mon <- vroom("dnmt3a_MON_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_bas <- vroom("dnmt3a_BAS_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_eos <- vroom("dnmt3a_EOS_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

dnmt3a_wbc <- rename_bcx_columns(dnmt3a_wbc, "wbc")
dnmt3a_neu <- rename_bcx_columns(dnmt3a_neu, "neu")
dnmt3a_lym <- rename_bcx_columns(dnmt3a_lym, "lym")
dnmt3a_mon <- rename_bcx_columns(dnmt3a_mon, "mon")
dnmt3a_bas <- rename_bcx_columns(dnmt3a_bas, "bas")
dnmt3a_eos <- rename_bcx_columns(dnmt3a_eos, "eos")

### Merge all datasets:
dnmt3a_bcx <- dnmt3a_wbc %>%
  inner_join(dnmt3a_neu, by = "SNP", suffix = c("_wbc", "_neu")) %>%
  inner_join(dnmt3a_lym, by = "SNP", suffix = c("", "_lym")) %>%
  inner_join(dnmt3a_mon, by = "SNP", suffix = c("", "_mon")) %>%
  inner_join(dnmt3a_bas, by = "SNP", suffix = c("", "_bas")) %>%
  inner_join(dnmt3a_eos, by = "SNP", suffix = c("", "_eos"))

### Remove duplicate rows:
dnmt3a_bcx <- dnmt3a_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- dnmt3a_bcx$SNP[duplicated(dnmt3a_bcx$SNP)]
duplicates

### Save output:
vroom_write(dnmt3a_bcx, file = "dnmt3a_bcx_wbc_harmonised_combined.tsv.gz")
