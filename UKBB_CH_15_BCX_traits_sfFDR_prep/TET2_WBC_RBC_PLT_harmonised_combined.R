### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
tet2_wbc <- vroom("tet2_bcx_wbc_harmonised_combined.tsv.gz")
tet2_rbc <- vroom("tet2_bcx_rbc_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2))
tet2_platelet <- vroom("tet2_bcx_platelet_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2))

### Combine:
tet2_bcx <- tet2_wbc %>%
  inner_join(tet2_rbc, by = "SNP") %>%
  inner_join(tet2_platelet, by = "SNP")
dim(tet2_bcx)

### Check for duplicates:
duplicates <- tet2_bcx$SNP[duplicated(tet2_bcx$SNP)]
duplicates

### Filter for p-value columns for sffdr:
tet2_bcx2 <- tet2_bcx %>%
  select("SNP", "rsid", "CHR", "BP", "pval_tet2", 
         "pval_wbc", "pval_neu", "pval_lym", "pval_mon","pval_bas", "pval_eos", 
         "pval_rbc", "pval_hgb", "pval_hct", "pval_mch", "pval_mcv", "pval_mchc", "pval_rdw", 
         "pval_plt", "pval_mpv")

 ### Save output:
vroom_write(tet2_bcx, file = "tet2_bcx_wbc_rbc_platelet_harmonised_combined.tsv.gz")
vroom_write(tet2_bcx2, file = "tet2_bcx_wbc_rbc_platelet_harmonised_combined_sffdr.tsv.gz")

### Separate harmonised stats by trait:
### Function to process summary statistics for a specific trait:
process_trait_data <- function(data, trait_suffix, N = 346787) {
  data %>%
    select(rsid, CHR, BP, 
           starts_with(paste0("ea_", trait_suffix)),
           starts_with(paste0("oa_", trait_suffix)),
           starts_with(paste0("eaf_", trait_suffix)),
           starts_with(paste0("beta_", trait_suffix)),
           starts_with(paste0("se_", trait_suffix)),
           starts_with(paste0("pval_", trait_suffix)),
           starts_with(paste0("n_samples_", trait_suffix))) %>%
    rename(SNP = rsid,
           A1 = starts_with(paste0("ea_", trait_suffix)),
           A2 = starts_with(paste0("oa_", trait_suffix)),
           EAF = starts_with(paste0("eaf_", trait_suffix)),
           B = starts_with(paste0("beta_", trait_suffix)),
           SE = starts_with(paste0("se_", trait_suffix)),
           P = starts_with(paste0("pval_", trait_suffix)),
           N = starts_with(paste0("n_samples_", trait_suffix)))}

### Process each blood cell trait:
traits <- c("tet2", "wbc", "neu", "lym", "mon", "bas", "eos", 
	        "rbc", "hgb", "hct", "mch", "mcv", "mchc", "rdw",
	        "plt", "mpv")
for (trait in traits) {
  processed_data <- process_trait_data(tet2_bcx, trait)
  output_file <- paste0(trait, "_tet2.sumstats.txt")
  write_delim(processed_data, file = output_file)}
