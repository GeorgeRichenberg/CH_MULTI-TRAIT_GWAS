### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_wbc <- vroom("overall_wbc_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_neu <- vroom("overall_neu_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_lym <- vroom("overall_lym_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_mon <- vroom("overall_mon_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_bas <- vroom("overall_bas_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_eos <- vroom("overall_eos_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

overall_wbc <- rename_bcx_columns(overall_wbc, "wbc")
overall_neu <- rename_bcx_columns(overall_neu, "neu")
overall_lym <- rename_bcx_columns(overall_lym, "lym")
overall_mon <- rename_bcx_columns(overall_mon, "mon")
overall_bas <- rename_bcx_columns(overall_bas, "bas")
overall_eos <- rename_bcx_columns(overall_eos, "eos")

### Merge all datasets:
overall_bcx <- overall_wbc %>%
  inner_join(overall_neu, by = "SNP", suffix = c("_wbc", "_neu")) %>%
  inner_join(overall_lym, by = "SNP", suffix = c("", "_lym")) %>%
  inner_join(overall_mon, by = "SNP", suffix = c("", "_mon")) %>%
  inner_join(overall_bas, by = "SNP", suffix = c("", "_bas")) %>%
  inner_join(overall_eos, by = "SNP", suffix = c("", "_eos"))

### Remove duplicate rows:
overall_bcx <- overall_bcx %>%
  distinct()

### Check for duplicates:  
duplicates <- overall_bcx$SNP[duplicated(overall_bcx$SNP)]
duplicates

### Save output:
vroom_write(overall_bcx, file = "overall_bcx_wbc_harmonised_combined.tsv.gz")