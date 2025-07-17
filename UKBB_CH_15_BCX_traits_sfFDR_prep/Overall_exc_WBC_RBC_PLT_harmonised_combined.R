### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_wbc <- vroom("overall_bcx_wbc_harmonised_combined.tsv.gz")
overall_rbc <- vroom("overall_bcx_rbc_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall))
overall_platelet <- vroom("overall_bcx_platelet_harmonised_combined.tsv.gz", col_select = -c(rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall))

### Combine:
overall_bcx <- overall_wbc %>%
  inner_join(overall_rbc, by = "SNP") %>%
  inner_join(overall_platelet, by = "SNP")
dim(overall_bcx)

### Check for duplicates:
duplicates <- overall_bcx$SNP[duplicated(overall_bcx$SNP)]
duplicates

### Filter for p-value columns for sffdr:
overall_bcx2 <- overall_bcx %>%
  select("SNP", "rsid", "CHR", "BP", "pval_overall", 
         "pval_wbc", "pval_neu", "pval_lym", "pval_mon","pval_bas", "pval_eos", 
         "pval_rbc", "pval_hgb", "pval_hct", "pval_mch", "pval_mcv", "pval_mchc", "pval_rdw", 
         "pval_plt", "pval_mpv")

 ### Save output:
vroom_write(overall_bcx, file = "overall_bcx_wbc_rbc_platelet_harmonised_combined.tsv.gz")
vroom_write(overall_bcx2, file = "overall_bcx_wbc_rbc_platelet_harmonised_combined_sffdr.tsv.gz")

### Separate harmonised stats by trait:
### Function to process summary statistics for a specific trait:
process_trait_data <- function(data, trait_suffix, N = 368526) {
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
traits <- c("overall", "wbc", "neu", "lym", "mon", "bas", "eos", 
	        "rbc", "hgb", "hct", "mch", "mcv", "mchc", "rdw",
	        "plt", "mpv")
for (trait in traits) {
  processed_data <- process_trait_data(overall_bcx, trait)
  output_file <- paste0(trait, "_overall.sumstats.txt")
  write_delim(processed_data, file = output_file)}
