### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_inc_sffdr_all <- vroom("overall_inc_ch_bcx_eur_sffdr_results.txt", col_select = c(SNP, Functional_pvalue))
overall_inc_kessler <- vroom("overall_GCST90165267_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall))

overall_exc_sffdr_all <- vroom("overall_ch_bcx_eur_sffdr_results.txt", col_select = c(SNP, Functional_pvalue))
overall_exc_kessler <- vroom("overall_GCST90165261_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall))

dnmt3a_sffdr_all <- vroom("dnmt3a_ch_bcx_eur_sffdr_results.txt", col_select = c(SNP, Functional_pvalue))
dnmt3a_kessler <- vroom("dnmt3a_GCST90165271_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a))

tet2_sffdr_all <- vroom("tet2_ch_bcx_eur_sffdr_results.txt", col_select = c(SNP, Functional_pvalue))
tet2_kessler <- vroom("tet2_GCST90165281_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2))

asxl1_sffdr_all <- vroom("asxl1_ch_bcx_eur_sffdr_results.txt", col_select = c(SNP, Functional_pvalue))
asxl1_kessler <- vroom("asxl1_GCST90165259_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(rsid, CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1))


### Prepare overall_inc sffdr fp for fuma:
overall_inc_kessler <- overall_inc_kessler %>%
  rename(SNP = rsid,
         A1 = ea_overall, 
         A2 = oa_overall,
         EAF = eaf_overall,
         Beta = beta_overall,
         SE = se_overall,
         P = pval_overall)
overall_inc_kessler_sffdr <- overall_inc_kessler %>%
  inner_join(overall_inc_sffdr_all, by = "SNP")  
overall_inc_kessler_sffdr <- overall_inc_kessler_sffdr %>%
  rename(FP = Functional_pvalue)
write_delim(overall_inc_kessler_sffdr, file = "overall_inc_kessler_sffdr_fp.txt.gz")

### Prepare overall_exc sffdr fp for fuma:
overall_exc_kessler <- overall_exc_kessler %>%
  rename(SNP = rsid,
         A1 = ea_overall, 
         A2 = oa_overall,
         EAF = eaf_overall,
         Beta = beta_overall,
         SE = se_overall,
         P = pval_overall)
overall_exc_kessler_sffdr <- overall_exc_kessler %>%
  inner_join(overall_exc_sffdr_all, by = "SNP")  
overall_exc_kessler_sffdr <- overall_exc_kessler_sffdr %>%
  rename(FP = Functional_pvalue)
write_delim(overall_exc_kessler_sffdr, file = "overall_exc_kessler_sffdr_fp.txt.gz")

### Prepare dnmt3a sffdr fp for fuma:
dnmt3a_kessler <- dnmt3a_kessler %>%
  rename(SNP = rsid,
         A1 = ea_dnmt3a, 
         A2 = oa_dnmt3a,
         EAF = eaf_dnmt3a,
         Beta = beta_dnmt3a,
         SE = se_dnmt3a,
         P = pval_dnmt3a)
dnmt3a_kessler_sffdr <- dnmt3a_kessler %>%
  inner_join(dnmt3a_sffdr_all, by = "SNP")  
dnmt3a_kessler_sffdr <- dnmt3a_kessler_sffdr %>%
  rename(FP = Functional_pvalue)
write_delim(dnmt3a_kessler_sffdr, file = "dnmt3a_kessler_sffdr_fp.txt.gz")

### Prepare tet2 sffdr fp for fuma:
tet2_kessler <- tet2_kessler %>%
  rename(SNP = rsid,
         A1 = ea_tet2, 
         A2 = oa_tet2,
         EAF = eaf_tet2,
         Beta = beta_tet2,
         SE = se_tet2,
         P = pval_tet2)
tet2_kessler_sffdr <- tet2_kessler %>%
  inner_join(tet2_sffdr_all, by = "SNP")  
tet2_kessler_sffdr <- tet2_kessler_sffdr %>%
  rename(FP = Functional_pvalue)
write_delim(tet2_kessler_sffdr, file = "tet2_kessler_sffdr_fp.txt.gz")

### Prepare asxl1 sffdr fp for fuma:
asxl1_kessler <- asxl1_kessler %>%
  rename(SNP = rsid,
         A1 = ea_asxl1, 
         A2 = oa_asxl1,
         EAF = eaf_asxl1,
         Beta = beta_asxl1,
         SE = se_asxl1,
         P = pval_asxl1)
asxl1_kessler_sffdr <- asxl1_kessler %>%
  inner_join(asxl1_sffdr_all, by = "SNP") 
asxl1_kessler_sffdr <- asxl1_kessler_sffdr %>%
  rename(FP = Functional_pvalue)
write_delim(asxl1_kessler_sffdr, file = "asxl1_kessler_sffdr_fp.txt.gz")