### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
tet2_rbc <- vroom("tet2_RBC_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_hgb <- vroom("tet2_HGB_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_hct <- vroom("tet2_HCT_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_mch <- vroom("tet2_MCH_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_mcv <- vroom("tet2_MCV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_mchc <- vroom("tet2_MCHC_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_rdw <- vroom("tet2_RDW_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

tet2_rbc <- rename_bcx_columns(tet2_rbc, "rbc")
tet2_hgb <- rename_bcx_columns(tet2_hgb, "hgb")
tet2_hct <- rename_bcx_columns(tet2_hct, "hct")
tet2_mch <- rename_bcx_columns(tet2_mch, "mch")
tet2_mcv <- rename_bcx_columns(tet2_mcv, "mcv")
tet2_mchc <- rename_bcx_columns(tet2_mchc, "mchc")
tet2_rdw <- rename_bcx_columns(tet2_rdw, "rdw")

### Merge all datasets:
tet2_bcx <- tet2_rbc %>%
  inner_join(tet2_hgb, by = "SNP", suffix = c("_rbc", "_hgb")) %>%
  inner_join(tet2_hct, by = "SNP", suffix = c("", "_hct")) %>%
  inner_join(tet2_mch, by = "SNP", suffix = c("", "_mch")) %>%
  inner_join(tet2_mcv, by = "SNP", suffix = c("", "_mcv")) %>%
  inner_join(tet2_mchc, by = "SNP", suffix = c("", "_mchc")) %>%
  inner_join(tet2_rdw, by = "SNP", suffix = c("", "_rdw"))

### Remove duplicate rows:
tet2_bcx <- tet2_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- tet2_bcx$SNP[duplicated(tet2_bcx$SNP)]
duplicates

### Save output:
vroom_write(tet2_bcx, file = "tet2_bcx_rbc_harmonised_combined.tsv.gz")