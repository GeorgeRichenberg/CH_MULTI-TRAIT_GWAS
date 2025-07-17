### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)

### Define a function to load, process, and save the dataset:
process_and_save <- function(input_file, output_file) {
  data <- vroom(input_file, col_select = c(rs_number, reference_allele, other_allele, eaf, beta, se, 'p-value', n_studies, n_samples))
    data <- data %>%
    filter(eaf >= 0.005, eaf <= 0.995) %>%
    filter(!reference_allele %in% c("I", "D"))
    data <- data %>%
    separate(rs_number, into = c("chr", "pos_ea_oa"), sep = ":", remove = FALSE) %>%
    separate(pos_ea_oa, into = c("bp", "ea_oa"), sep = "_", remove = TRUE) %>%
    select(-ea_oa)
    data <- data %>%
    rename(ea_bcx = reference_allele,
           oa_bcx = other_allele,
           eaf_bcx = eaf,
           beta_bcx = beta,
           se_bcx = se,
           pval_bcx = `p-value`,
           n_studies_bcx = n_studies,
           n_samples_bcx = n_samples) %>%
    distinct()
    vroom_write(data, file = output_file)}

### Define input and output file names:
files <- list(
  list(input = "RBC_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "RBC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "HGB_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "HGB_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "HCT_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "HCT_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "MCH_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "MCH_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "MCV_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "MCV_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "MCHC_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "MCHC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "RDW_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "RDW_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),)

### Loop through each file and apply the process_and_save function:
for (file in files) {
  process_and_save(file$input, file$output)}