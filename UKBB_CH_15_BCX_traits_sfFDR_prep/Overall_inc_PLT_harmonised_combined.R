### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_plt <- vroom("overall_inc_PLT_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
overall_mpv <- vroom("overall_inc_MPV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

overall_plt <- rename_bcx_columns(overall_plt, "plt")
overall_mpv <- rename_bcx_columns(overall_mpv, "mpv")

### Merge the datasets on SNP:
overall_bcx <- overall_plt %>%
  inner_join(overall_mpv, by = "SNP", suffix = c("_plt", "_mpv"))

### Remove duplicate rows:
overall_bcx <- overall_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- overall_bcx$SNP[duplicated(overall_bcx$SNP)]
duplicates

### Save output:
vroom_write(overall_bcx, file = "overall_inc_bcx_platelet_harmonised_combined.tsv.gz")