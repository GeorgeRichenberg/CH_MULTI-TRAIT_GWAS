### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

# Load data:
overall <- vroom("overall_GCST90165267_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
blood_cell_traits <- list(
  WBC = "WBC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  NEU = "NEU_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  LYM = "LYM_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  MON = "MON_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  BAS = "BAS_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  EOS = "EOS_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz")

### Create function to process and save harmonised and unharmonised outputs for each trait:
process_trait <- function(trait, input_file) {
  trait_data <- vroom(input_file)
  overall <- overall %>%
    filter(str_length(ea_overall) == 1 & str_length(oa_overall) == 1)
  trait_data <- trait_data %>%
    mutate(SNP = paste(chr, bp, sep = "_"))
  merged_data <- overall %>%
    inner_join(trait_data, by = "SNP")
  mismatched_alleles <- merged_data %>%
    rowwise() %>%
    filter(!(ea_overall %in% c(ea_bcx, oa_bcx) & oa_overall %in% c(ea_bcx, oa_bcx))) %>%
    ungroup()
  merged_data <- merged_data %>%
    anti_join(mismatched_alleles, by = c("SNP", "ea_overall", "oa_overall", "ea_bcx", "oa_bcx"))
  merged_data <- merged_data %>%
    group_by(SNP) %>%
    filter(pval_overall == min(pval_overall)) %>%
    ungroup()
  harmonised_data <- merged_data %>%
    mutate(new_beta_bcx = ifelse(ea_overall != ea_bcx & oa_overall != oa_bcx,
                             beta_bcx * -1,
                             beta_bcx),
           new_eaf_bcx = ifelse(ea_overall != ea_bcx & oa_overall != oa_bcx,
                            1 - eaf_bcx,
                            eaf_bcx),
           temp_ea_bcx = ea_bcx,
           temp_oa_bcx = oa_bcx,
           new_ea_bcx = ifelse(ea_overall != ea_bcx & oa_overall != oa_bcx, 
                           temp_oa_bcx, 
                           ea_bcx),
           new_oa_bcx = ifelse(ea_overall != ea_bcx & oa_overall != oa_bcx, 
                           temp_ea_bcx, 
                           oa_bcx)) %>%
    select(-temp_ea_bcx, -temp_oa_bcx)
  harmonised_data <- harmonised_data %>%
    select(-c(rs_number, chr, bp)) %>%
    select(SNP, rsid, CHR, BP,
           ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall, 
           ea_bcx, oa_bcx, new_ea_bcx, new_oa_bcx, eaf_bcx, new_eaf_bcx, beta_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx)
  harmonised_file <- paste0("overall_inc_", trait, "_harmonised.tsv.gz")
  vroom_write(harmonised_data, file = harmonised_file)}

### Apply the processing function to each trait:
for (trait in names(blood_cell_traits)) {
  process_trait(trait, blood_cell_traits[[trait]])}
