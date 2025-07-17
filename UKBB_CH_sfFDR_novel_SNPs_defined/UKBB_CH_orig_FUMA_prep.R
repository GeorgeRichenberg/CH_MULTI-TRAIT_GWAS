### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_inc_kessler <- vroom("overall_GCST90165267_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_overall, oa_overall, beta_overall, se_overall, pval_overall))
overall_exc_kessler <- vroom("overall_GCST90165261_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_overall, oa_overall, beta_overall, se_overall, pval_overall))
dnmt3a_kessler <- vroom("dnmt3a_GCST90165271_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a))
tet2_kessler <- vroom("tet2_GCST90165281_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_tet2, oa_tet2, beta_tet2, se_tet2, pval_tet2))
asxl1_kessler <- vroom("asxl1_GCST90165259_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_asxl1, oa_asxl1, beta_asxl1, se_asxl1, pval_asxl1))

### Prepare overall_inc kessler for fuma:
overall_inc_kessler <- overall_inc_kessler %>%
  rename(SNP = rsid,
         A1 = ea_overall, 
         A2 = oa_overall,
         Beta = beta_overall,
         SE = se_overall,
         P = pval_overall)
write_delim(overall_inc_kessler, file = "overall_inc_kessler_orig_fuma.txt.gz")

### Prepare overall_exc kessler for fuma:
overall_exc_kessler <- overall_exc_kessler %>%
  rename(SNP = rsid,
         A1 = ea_overall, 
         A2 = oa_overall,
         Beta = beta_overall,
         SE = se_overall,
         P = pval_overall)
write_delim(overall_exc_kessler, file = "overall_exc_kessler_orig_fuma.txt.gz")

### Prepare dnmt3a kessler for fuma:
dnmt3a_kessler <- dnmt3a_kessler %>%
  rename(SNP = rsid,
         A1 = ea_dnmt3a, 
         A2 = oa_dnmt3a,
         Beta = beta_dnmt3a,
         SE = se_dnmt3a,
         P = pval_dnmt3a)
write_delim(dnmt3a_kessler, file = "dnmt3a_kessler_orig_fuma.txt.gz")

### Prepare tet2 kessler for fuma:
tet2_kessler <- tet2_kessler %>%
  rename(SNP = rsid,
         A1 = ea_tet2, 
         A2 = oa_tet2,
         Beta = beta_tet2,
         SE = se_tet2,
         P = pval_tet2)
write_delim(tet2_kessler, file = "tet2_kessler_orig_fuma.txt.gz")

### Prepare asxl1 kessler for fuma:
asxl1_kessler <- asxl1_kessler %>%
  rename(SNP = rsid,
         A1 = ea_asxl1, 
         A2 = oa_asxl1,
         Beta = beta_asxl1,
         SE = se_asxl1,
         P = pval_asxl1)
write_delim(asxl1_kessler, file = "asxl1_kessler_orig_fuma.txt.gz")