### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, p, fp))
mcaaut <- vroom("mCAaut_GCST90165285_buildGRCh37_selectedcols_1pcMAF_rsid_beta_se.tsv.gz")
mlox <- vroom("mLOX_GCST90165287_buildGRCh37_selectedcols_1pcMAF_rsid_beta_se.tsv.gz")
mloy <- vroom("mLOY_GCST90165289_buildGRCh37_selectedcols_1pcMAF_rsid_beta_se.tsv.gz")

### Format sffdr_novel:
sffdr_novel <- sffdr_novel %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR|\\[|\\]"))
sffdr_novel2 <- sffdr_novel %>%
  mutate(id = paste0(rsID, "_", pmin(ea, oa), "_", pmax(ea, oa))) %>%
  select(id, everything())

### Format mca:
mca <- mcaaut %>%
  mutate(id = paste0(rsid , "_", pmin(effect_allele, other_allele), "_", pmax(effect_allele, other_allele))) %>%
  rename(ea_mca = effect_allele,
         oa_mca = other_allele,
         eaf_mca = effect_allele_frequency,
         beta_mca = beta,
         se_mca = se, 
         pval_mca = p_value) %>%
  select(id, ea_mca, oa_mca, eaf_mca, beta_mca, se_mca, pval_mca)
mlox <- mlox %>%
  mutate(id = paste0(rsid , "_", pmin(effect_allele, other_allele), "_", pmax(effect_allele, other_allele))) %>%
  rename(ea_mlox = effect_allele,
         oa_mlox = other_allele,
         eaf_mlox = effect_allele_frequency,
         beta_mlox = beta,
         se_mlox = se, 
         pval_mlox = p_value) %>%
  select(id, ea_mlox, oa_mlox, eaf_mlox, beta_mlox, se_mlox, pval_mlox)
mloy <- mloy %>%
  mutate(id = paste0(rsid , "_", pmin(effect_allele, other_allele), "_", pmax(effect_allele, other_allele))) %>%
  rename(ea_mloy = effect_allele,
         oa_mloy = other_allele,
         eaf_mloy = effect_allele_frequency,
         beta_mloy = beta,
         se_mloy = se, 
         pval_mloy = p_value) %>%
  select(id, ea_mloy, oa_mloy, eaf_mloy, beta_mloy, se_mloy, pval_mloy)

### Merge sffdr_novel with mosaic data:
sffdr_mca <- sffdr_novel2 %>%
  inner_join(mca, by = "id")
dim(sffdr_mca)
sffdr_mlox <- sffdr_novel2 %>%
  inner_join(mlox, by = "id")
dim(sffdr_mlox)
sffdr_mloy <- sffdr_novel2 %>%
  inner_join(mloy, by = "id")
dim(sffdr_mloy)
sffdr_mos <- sffdr_mca %>%
  inner_join(sffdr_mlox, by = c("id", "CH_GWAS", "rsID", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p", "fp")) %>%
  inner_join(sffdr_mloy, by = c("id", "CH_GWAS", "rsID", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p", "fp"))
dim(sffdr_mos)

### Harmonise data to CH risk increasing allele:
sffdr_mos <- sffdr_mos %>%
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
    
    beta_mca = ifelse(flip, -beta_mca, beta_mca),
    eaf_mca = ifelse(flip, 1 - eaf_mca, eaf_mca),
    or_mca = exp(beta_mca),
    l95cl_mca = exp(beta_mca - 1.96 * se_mca),
    u95cl_mca = exp(beta_mca + 1.96 * se_mca),
    
    beta_mlox = ifelse(flip, -beta_mlox, beta_mlox),
    eaf_mlox = ifelse(flip, 1 - eaf_mlox, eaf_mlox),
    or_mlox = exp(beta_mlox),
    l95cl_mlox = exp(beta_mlox - 1.96 * se_mlox),
    u95cl_mlox = exp(beta_mlox + 1.96 * se_mlox),
    
    beta_mloy = ifelse(flip, -beta_mloy, beta_mloy),
    eaf_mloy = ifelse(flip, 1 - eaf_mloy, eaf_mloy),
    or_mloy = exp(beta_mloy),
    l95cl_mloy = exp(beta_mloy - 1.96 * se_mloy),
    u95cl_mloy = exp(beta_mloy + 1.96 * se_mloy)) %>%
  select(-flip, -ea_orig, -oa_orig)

sffdr_mos  <- sffdr_mos  %>%
  select(id, CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp,
         eaf_mca, beta_mca, se_mca, or_mca, l95cl_mca, u95cl_mca, pval_mca,
         eaf_mlox, beta_mlox, se_mlox, or_mlox, l95cl_mlox, u95cl_mlox, pval_mlox,
         eaf_mloy, beta_mloy, se_mloy, or_mloy, l95cl_mloy, u95cl_mloy, pval_mloy)

### Save:
write_csv(sffdr_mos, file = "sffdr_novel_mca_mlox_mloy.csv")