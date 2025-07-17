### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
asxl1_plt <- vroom("asxl1_PLT_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
asxl1_mpv <- vroom("asxl1_MPV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

asxl1_plt <- rename_bcx_columns(asxl1_plt, "plt")
asxl1_mpv <- rename_bcx_columns(asxl1_mpv, "mpv")

### Merge the datasets on SNP:
asxl1_bcx <- asxl1_plt %>%
  inner_join(asxl1_mpv, by = "SNP", suffix = c("_plt", "_mpv"))

### Remove duplicate rows:
asxl1_bcx <- asxl1_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- asxl1_bcx$SNP[duplicated(asxl1_bcx$SNP)]
duplicates

### Save output:
vroom_write(asxl1_bcx, file = "asxl1_bcx_platelet_harmonised_combined.tsv.gz")