### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
tet2_plt <- vroom("tet2_PLT_harmonised.tsv.gz", col_select = c(SNP, rsid, CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))
tet2_mpv <- vroom("tet2_MPV_harmonised.tsv.gz", col_select = c(SNP, new_ea_bcx, new_oa_bcx, new_eaf_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx))

### Rename columns with respective blood cell traits:
rename_bcx_columns <- function(data, trait_suffix) {
  colnames(data) <- gsub("^new_(.*)_bcx$", paste0("\\1_", trait_suffix), colnames(data))
  colnames(data) <- gsub("_bcx$", paste0("_", trait_suffix), colnames(data))
  return(data)}

tet2_plt <- rename_bcx_columns(tet2_plt, "plt")
tet2_mpv <- rename_bcx_columns(tet2_mpv, "mpv")

### Merge the datasets on SNP:
tet2_bcx <- tet2_plt %>%
  inner_join(tet2_mpv, by = "SNP", suffix = c("_plt", "_mpv"))

### Remove duplicate rows:
tet2_bcx <- tet2_bcx %>%
  distinct()

### Check for duplicates:
duplicates <- tet2_bcx$SNP[duplicated(tet2_bcx$SNP)]
duplicates

### Save output:
vroom_write(tet2_bcx, file = "tet2_bcx_platelet_harmonised_combined.tsv.gz")