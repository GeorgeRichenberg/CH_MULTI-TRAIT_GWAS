### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
tet2_wbc <- vroom("tet2_wbc_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_neu <- vroom("tet2_neu_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_lym <- vroom("tet2_lym_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_mon <- vroom("tet2_mon_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_bas <- vroom("tet2_bas_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_eos <- vroom("tet2_eos_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

tet2_wbc <- rename_bcx_columns(tet2_wbc, "wbc")
tet2_neu <- rename_bcx_columns(tet2_neu, "neu")
tet2_lym <- rename_bcx_columns(tet2_lym, "lym")
tet2_mon <- rename_bcx_columns(tet2_mon, "mon")
tet2_bas <- rename_bcx_columns(tet2_bas, "bas")
tet2_eos <- rename_bcx_columns(tet2_eos, "eos")

### Merge all datasets:
tet2_bcx <- tet2_wbc %>%
  inner_join(tet2_neu, by = "SNP", suffix = c("_wbc", "_neu")) %>%
  inner_join(tet2_lym, by = "SNP", suffix = c("", "_lym")) %>%
  inner_join(tet2_mon, by = "SNP", suffix = c("", "_mon")) %>%
  inner_join(tet2_bas, by = "SNP", suffix = c("", "_bas")) %>%
  inner_join(tet2_eos, by = "SNP", suffix = c("", "_eos"))

### Remove duplicate rows:
tet2_bcx <- tet2_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- tet2_bcx$SNP[duplicated(tet2_bcx$SNP)]
duplicates

### Save output:
vroom_write(tet2_bcx, file = "tet2_bcx_wbc_harmonised_combined.tsv.gz")