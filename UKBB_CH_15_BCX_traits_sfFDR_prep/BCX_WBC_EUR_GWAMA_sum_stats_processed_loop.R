### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)

# Define a function to load, process, and save the dataset
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

### Define the input and output files:
files <- list(
  list(input = "WBC_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "WBC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "NEU_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "NEU_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "LYM_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "LYM_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "MON_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "MON_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "BAS_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "BAS_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
  list(input = "EOS_Europeans_GWAMA_0205_2019_newINTERVAL_out.out", output = "EOS_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"))

### Loop through each file and apply the process_and_save function:
for (file in files) {
  process_and_save(file$input, file$output)}
