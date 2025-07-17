### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
dnmt3a_wbc <- vroom("dnmt3a_bcx_wbc_harmonised_combined.tsv.gz")
dnmt3a_rbc <- vroom("dnmt3a_bcx_rbc_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a))
dnmt3a_platelet <- vroom("dnmt3a_bcx_platelet_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a))

### Combine:
dnmt3a_bcx <- dnmt3a_wbc %>%
  inner_join(dnmt3a_rbc, by = "SNP") %>%
  inner_join(dnmt3a_platelet, by = "SNP")
dim(dnmt3a_bcx)

### Check for duplicates:
duplicates <- dnmt3a_bcx$SNP[duplicated(dnmt3a_bcx$SNP)]
duplicates

### Filter for p-value columns for sffdr:
dnmt3a_bcx2 <- dnmt3a_bcx %>%
  select("SNP", "rsid", "CHR", "BP", "pval_dnmt3a", 
         "pval_wbc", "pval_neu", "pval_lym", "pval_mon","pval_bas", "pval_eos", 
         "pval_rbc", "pval_hgb", "pval_hct", "pval_mch", "pval_mcv", "pval_mchc", "pval_rdw", 
         "pval_plt", "pval_mpv")

 ### Save output:
vroom_write(dnmt3a_bcx, file = "dnmt3a_bcx_wbc_rbc_platelet_harmonised_combined.tsv.gz")
vroom_write(dnmt3a_bcx2, file = "dnmt3a_bcx_wbc_rbc_platelet_harmonised_combined_sffdr.tsv.gz")

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
traits <- c("dnmt3a", "wbc", "neu", "lym", "mon", "bas", "eos", 
	        "rbc", "hgb", "hct", "mch", "mcv", "mchc", "rdw",
	        "plt", "mpv")
for (trait in traits) {
  processed_data <- process_trait_data(dnmt3a_bcx, trait)
  output_file <- paste0(trait, "_dnmt3a.sumstats.txt")
  write_delim(processed_data, file = output_file)}
