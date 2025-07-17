### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
asxl1_wbc <- vroom("asxl1_wbc_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_neu <- vroom("asxl1_neu_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_lym <- vroom("asxl1_lym_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_mon <- vroom("asxl1_mon_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_bas <- vroom("asxl1_bas_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_eos <- vroom("asxl1_eos_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

asxl1_wbc <- rename_bcx_columns(asxl1_wbc, "wbc")
asxl1_neu <- rename_bcx_columns(asxl1_neu, "neu")
asxl1_lym <- rename_bcx_columns(asxl1_lym, "lym")
asxl1_mon <- rename_bcx_columns(asxl1_mon, "mon")
asxl1_bas <- rename_bcx_columns(asxl1_bas, "bas")
asxl1_eos <- rename_bcx_columns(asxl1_eos, "eos")

### Merge all datasets:
asxl1_bcx <- asxl1_wbc %>%
  inner_join(asxl1_neu, by = "SNP", suffix = c("_wbc", "_neu")) %>%
  inner_join(asxl1_lym, by = "SNP", suffix = c("", "_lym")) %>%
  inner_join(asxl1_mon, by = "SNP", suffix = c("", "_mon")) %>%
  inner_join(asxl1_bas, by = "SNP", suffix = c("", "_bas")) %>%
  inner_join(asxl1_eos, by = "SNP", suffix = c("", "_eos"))

### Remove duplicate rows:
asxl1_bcx <- asxl1_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- asxl1_bcx$SNP[duplicated(asxl1_bcx$SNP)]
duplicates

### Save output:
vroom_write(asxl1_bcx, file = "asxl1_bcx_wbc_harmonised_combined.tsv.gz")