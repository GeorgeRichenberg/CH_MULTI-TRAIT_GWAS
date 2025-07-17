### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, p, fp))
ltl_pc1 <- vroom("GCST90435144.tsv.gz", col_select = c(rs_id, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value))
ltl_pc2 <- vroom("GCST90435145.tsv.gz", col_select = c(rs_id, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value))

### Format sffdr_novel:
sffdr_novel <- sffdr_novel %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR|\\[|\\]"))
sffdr_novel2 <- sffdr_novel %>%
  mutate(id = paste0(rsID, "_", pmin(ea, oa), "_", pmax(ea, oa))) %>%
  select(id, everything())

### Format ltl:
ltl_pc1 <- ltl_pc1 %>%
  mutate(id = paste0(rs_id , "_", pmin(effect_allele, other_allele), "_", pmax(effect_allele, other_allele))) %>%
  rename(ea_ltl_pc1 = effect_allele,
         oa_ltl_pc1 = other_allele,
         eaf_ltl_pc1 = effect_allele_frequency,
         beta_ltl_pc1 = beta,
         se_ltl_pc1 = standard_error, 
         pval_ltl_pc1 = p_value) %>%
  select(id, ea_ltl_pc1, oa_ltl_pc1, eaf_ltl_pc1, beta_ltl_pc1, se_ltl_pc1, pval_ltl_pc1)

ltl_pc2 <- ltl_pc2 %>%
  mutate(id = paste0(rs_id , "_", pmin(effect_allele, other_allele), "_", pmax(effect_allele, other_allele))) %>%
  rename(ea_ltl_pc2 = effect_allele,
         oa_ltl_pc2 = other_allele,
         eaf_ltl_pc2 = effect_allele_frequency,
         beta_ltl_pc2 = beta,
         se_ltl_pc2 = standard_error, 
         pval_ltl_pc2 = p_value) %>%
  select(id, ea_ltl_pc2, oa_ltl_pc2, eaf_ltl_pc2, beta_ltl_pc2, se_ltl_pc2, pval_ltl_pc2)

### Merge sffdr_novel with ltl data:
sffdr_ltl_pc1 <- sffdr_novel2 %>%
  left_join(ltl_pc1, by = "id")
dim(sffdr_ltl_pc1)
sffdr_ltl_pc2 <- sffdr_novel2 %>%
  left_join(ltl_pc2, by = "id")
dim(sffdr_ltl_pc2)
sffdr_ltl_pc <- sffdr_ltl_pc1 %>%
  inner_join(sffdr_ltl_pc2, by = c("id", "CH_GWAS", "rsID", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p", "fp"))

### Harmonise data to CH risk increasing allele:
sffdr_ltl_pc <- sffdr_ltl_pc %>%
  mutate(
    flip = beta < 0,
    ea_orig = ea,
    oa_orig = oa,
    
    ea = ifelse(flip, oa_orig, ea),
    oa = ifelse(flip, ea_orig, oa),
    beta = ifelse(flip, -beta, beta),
    eaf = ifelse(flip, 1 - eaf, eaf),
    or = exp(beta),
    l95cl = exp(beta - 1.96 * se),
    u95cl = exp(beta + 1.96 * se),
    
    beta_ltl_pc1 = ifelse(flip, -beta_ltl_pc1, beta_ltl_pc1),
    eaf_ltl_pc1 = ifelse(flip, 1 - eaf_ltl_pc1, eaf_ltl_pc1),
    or_ltl_pc1 = exp(beta_ltl_pc1),
    l95cl_ltl_pc1 = exp(beta_ltl_pc1 - 1.96 * se_ltl_pc1),
    u95cl_ltl_pc1 = exp(beta_ltl_pc1 + 1.96 * se_ltl_pc1),
    
    beta_ltl_pc2 = ifelse(flip, -beta_ltl_pc2, beta_ltl_pc2),
    eaf_ltl_pc2 = ifelse(flip, 1 - eaf_ltl_pc2, eaf_ltl_pc2),
    or_ltl_pc2 = exp(beta_ltl_pc2),
    l95cl_ltl_pc2 = exp(beta_ltl_pc2 - 1.96 * se_ltl_pc2),
    u95cl_ltl_pc2 = exp(beta_ltl_pc2 + 1.96 * se_ltl_pc2)) %>%
  select(-flip, -ea_orig, -oa_orig)

sffdr_ltl_pc  <- sffdr_ltl_pc  %>%
  select(id, CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp,
         eaf_ltl_pc1, beta_ltl_pc1, se_ltl_pc1, or_ltl_pc1, l95cl_ltl_pc1, u95cl_ltl_pc1, pval_ltl_pc1,
         eaf_ltl_pc2, beta_ltl_pc2, se_ltl_pc2, or_ltl_pc2, l95cl_ltl_pc2, u95cl_ltl_pc2, pval_ltl_pc2)

### Save:
write_csv(sffdr_ltl_pc, file = "sffdr_novel_ltl.csv")