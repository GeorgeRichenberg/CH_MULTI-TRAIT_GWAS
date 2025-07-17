### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
asxl1_wbc <- vroom("asxl1_bcx_wbc_harmonised_combined.tsv.gz")
asxl1_rbc <- vroom("asxl1_bcx_rbc_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1))
asxl1_platelet <- vroom("asxl1_bcx_platelet_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1))

### Combine:
asxl1_bcx <- asxl1_wbc %>%
  inner_join(asxl1_rbc, by = "SNP") %>%
  inner_join(asxl1_platelet, by = "SNP")
dim(asxl1_bcx)

### Check for duplicates:
duplicates <- asxl1_bcx$SNP[duplicated(asxl1_bcx$SNP)]
duplicates

### Filter for p-value columns for sffdr:
asxl1_bcx2 <- asxl1_bcx %>%
  select("SNP", "rsid", "CHR", "BP", "pval_asxl1", 
         "pval_wbc", "pval_neu", "pval_lym", "pval_mon","pval_bas", "pval_eos", 
         "pval_rbc", "pval_hgb", "pval_hct", "pval_mch", "pval_mcv", "pval_mchc", "pval_rdw", 
         "pval_plt", "pval_mpv")

 ### Save output:
vroom_write(asxl1_bcx, file = "asxl1_bcx_wbc_rbc_platelet_harmonised_combined.tsv.gz")
vroom_write(asxl1_bcx2, file = "asxl1_bcx_wbc_rbc_platelet_harmonised_combined_sffdr.tsv.gz")

### Separate harmonised stats by trait:
### Function to process summary statistics for a specific trait:
process_trait_data <- function(data, trait_suffix, N = 359088) {
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
traits <- c("asxl1", "wbc", "neu", "lym", "mon", "bas", "eos", 
	        "rbc", "hgb", "hct", "mch", "mcv", "mchc", "rdw",
	        "plt", "mpv")
for (trait in traits) {
  processed_data <- process_trait_data(asxl1_bcx, trait)
  output_file <- paste0(trait, "_asxl1.sumstats.txt")
  write_delim(processed_data, file = output_file)}
