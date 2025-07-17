### Set working directory and load required packages
setwd("")
library(vroom)
library(tidyverse)

### Load data:
asxl1 <- vroom("asxl1_GCST90165259_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
blood_cell_traits <- list(
  RBC = "RBC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  HGB = "HGB_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  HCT = "HCT_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  MCH = "MCH_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  MCV = "MCV_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  MCHC = "MCHC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz",
  RDW = "RDW_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz")

### Create function to process and save harmonised and unharmonised outputs for each trait:
process_trait <- function(trait, input_file) {
  trait_data <- vroom(input_file)
  asxl1_filtered <- asxl1 %>%
    filter(str_length(ea_asxl1) == 1 & str_length(oa_asxl1) == 1)
  trait_data <- trait_data %>%
    mutate(SNP = paste(chr, bp, sep = "_"))
  merged_data <- asxl1_filtered %>%
    inner_join(trait_data, by = "SNP")
  mismatched_alleles <- merged_data %>%
    rowwise() %>%
    filter(!(ea_asxl1 %in% c(ea_bcx, oa_bcx) & oa_asxl1 %in% c(ea_bcx, oa_bcx))) %>%
    ungroup()
  merged_data <- merged_data %>%
    anti_join(mismatched_alleles, by = c("SNP", "ea_asxl1", "oa_asxl1", "ea_bcx", "oa_bcx"))
  merged_data <- merged_data %>%
    group_by(SNP) %>%
    filter(pval_asxl1 == min(pval_asxl1)) %>%
    ungroup()
  harmonised_data <- merged_data %>%
    mutate(new_beta_bcx = ifelse(ea_asxl1 != ea_bcx & oa_asxl1 != oa_bcx,
                             beta_bcx * -1,
                             beta_bcx),
           new_eaf_bcx = ifelse(ea_asxl1 != ea_bcx & oa_asxl1 != oa_bcx,
                            1 - eaf_bcx,
                            eaf_bcx),
           temp_ea_bcx = ea_bcx,
           temp_oa_bcx = oa_bcx,
           new_ea_bcx = ifelse(ea_asxl1 != ea_bcx & oa_asxl1 != oa_bcx, 
                           temp_oa_bcx, 
                           ea_bcx),
           new_oa_bcx = ifelse(ea_asxl1 != ea_bcx & oa_asxl1 != oa_bcx, 
                           temp_ea_bcx, 
                           oa_bcx)) %>%
    select(-temp_ea_bcx, -temp_oa_bcx)
  harmonised_data <- harmonised_data %>%
    select(-c(rs_number, chr, bp)) %>%
    select(SNP, rsid, CHR, BP,
           ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1, 
           ea_bcx, oa_bcx, new_ea_bcx, new_oa_bcx, eaf_bcx, new_eaf_bcx, beta_bcx, new_beta_bcx, se_bcx, pval_bcx, n_studies_bcx, n_samples_bcx)
  harmonised_file <- paste0("asxl1_", trait, "_harmonised.tsv.gz")
  vroom_write(harmonised_data, file = harmonised_file)}

### Apply the processing function to each trait:
for (trait in names(blood_cell_traits)) {
  process_trait(trait, blood_cell_traits[[trait]])}