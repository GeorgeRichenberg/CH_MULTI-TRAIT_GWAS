### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
dnmt3a_plt <- vroom("dnmt3a_PLT_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
dnmt3a_mpv <- vroom("dnmt3a_MPV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}
  
dnmt3a_plt <- rename_bcx_columns(dnmt3a_plt, "plt")
dnmt3a_mpv <- rename_bcx_columns(dnmt3a_mpv, "mpv")

### Merge the datasets on SNP:
dnmt3a_bcx <- dnmt3a_plt %>%
  inner_join(dnmt3a_mpv, by = "SNP", suffix = c("_plt", "_mpv"))

### Remove duplicate rows:
dnmt3a_bcx <- dnmt3a_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- dnmt3a_bcx$SNP[duplicated(dnmt3a_bcx$SNP)]
duplicates

### Save output:
vroom_write(dnmt3a_bcx, file = "dnmt3a_bcx_platelet_harmonised_combined.tsv.gz")
