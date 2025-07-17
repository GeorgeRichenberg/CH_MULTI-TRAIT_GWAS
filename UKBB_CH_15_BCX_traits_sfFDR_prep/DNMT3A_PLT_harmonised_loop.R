### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
dnmt3a <- vroom("dnmt3a_GCST90165271_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
blood_cell_traits <- list(
  PLT = "PLT_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  MPV = "MPV_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz")


### Create function to process and save harmonised and unharmonised outputs for each trait:
process_trait <- function(trait, input_file) {
  trait_data <- vroom(input_file)
  dnmt3a_filtered <- dnmt3a %>%
    filter(str_length(ea_dnmt3a) == 1 & str_length(oa_dnmt3a) == 1)
  trait_data <- trait_data %>%
    mutate(SNP = paste(chr, bp, sep = "_"))
  merged_data <- dnmt3a_filtered %>%
    inner_join(trait_data, by = "SNP")
  mismatched_alleles <- merged_data %>%
    rowwise() %>%
    filter(!(ea_dnmt3a %in% c(ea_bcx, oa_bcx) & oa_dnmt3a %in% c(ea_bcx, oa_bcx))) %>%
    ungroup()
  merged_data <- merged_data %>%
    anti_join(mismatched_alleles, by = c("SNP", "ea_dnmt3a", "oa_dnmt3a", "ea_bcx", "oa_bcx"))
  merged_data <- merged_data %>%
    group_by(SNP) %>%
    filter(pval_dnmt3a == min(pval_dnmt3a)) %>%
    ungroup()
  harmonised_data <- merged_data %>%
    mutate(new_beta_bcx = ifelse(ea_dnmt3a != ea_bcx & oa_dnmt3a != oa_bcx,
                             beta_bcx * -1,
                             beta_bcx),
           new_eaf_bcx = ifelse(ea_dnmt3a != ea_bcx & oa_dnmt3a != oa_bcx,
                            1 - eaf_bcx,
                            eaf_bcx),
           temp_ea_bcx = ea_bcx,
           temp_oa_bcx = oa_bcx,
           new_ea_bcx = ifelse(ea_dnmt3a != ea_bcx & oa_dnmt3a != oa_bcx, 
                           temp_oa_bcx, 
                           ea_bcx),
           new_oa_bcx = ifelse(ea_dnmt3a != ea_bcx & oa_dnmt3a != oa_bcx, 
                           temp_ea_bcx, 
                           oa_bcx)) %>%
    select(-temp_ea_bcx, -temp_oa_bcx)
  harmonised_data <- harmonised_data %>%
    select(-c(rs_number, chr, bp)) %>%
    select(SNP, rsid, CHR, BP,
           ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a, 
           ea_bcx, oa_bcx, new_ea_bcx, new_oa_bcx, eaf_bcx, new_eaf_bcx, beta_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx)
  harmonised_file <- paste0("dnmt3a_", trait, "_harmonised.tsv.gz")
  vroom_write(harmonised_data, file = harmonised_file)}

### Apply the processing function to each trait:
for (trait in names(blood_cell_traits)) {
  process_trait(trait, blood_cell_traits[[trait]])}