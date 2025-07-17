### Set wd and load required packages:
setwd("")
library(vroom)
library(tidyverse)

#### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, p, fp))
aml <- vroom("finngen_R12_C3_AML_EXALLC.gz")
myelodys <- vroom("finngen_R12_C3_MYELODYSP_SYNDR_EXALLC.gz")
myeloprof <- vroom("finngen_R12_MYELOPROF_NONCML.gz")
myeloid <- vroom("finngen_R12_CD2_MYELOID_LEUKAEMIA_EXALLC.gz")

### Format sffdr_novel:
sffdr_novel <- sffdr_novel %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR|\\[|\\]"))
sffdr_novel2 <- sffdr_novel %>%
  mutate(id = paste0(rsID, "_", pmin(ea, oa), "_", pmax(ea, oa))) %>%
  select(id, everything())
rsid <- sffdr_novel2 %>%
  select(id)

### Format blood diseases:
aml <- aml %>%
  mutate(id = paste0(rsids, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  select(id, everything()) %>%
  inner_join(rsid, by = "id")
dim(aml)
aml <- aml %>%
  select(id, rsids, alt, ref, af_alt, beta, sebeta, pval) %>%
  rename(ea_aml = alt,
         oa_aml = ref,
         eaf_aml = af_alt,
         beta_aml = beta,
         se_aml = sebeta,
         p_aml = pval) %>%
  select(id, rsids, ea_aml, oa_aml, eaf_aml, beta_aml, se_aml, p_aml)

myelodys <- myelodys %>%
  mutate(id = paste0(rsids, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  select(id, everything()) %>%
  inner_join(rsid, by = "id")
dim(myelodys)
myelodys <- myelodys %>%
  select(id, rsids, alt, ref, af_alt, beta, sebeta, pval) %>%
  rename(ea_myelodys = alt,
         oa_myelodys = ref,
         eaf_myelodys = af_alt,
         beta_myelodys = beta,
         se_myelodys = sebeta,
         p_myelodys = pval) %>%
  select(id, rsids, ea_myelodys, oa_myelodys, eaf_myelodys, beta_myelodys, se_myelodys, p_myelodys)

myeloprof <- myeloprof %>%
  mutate(id = paste0(rsids, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  select(id, everything()) %>%
  inner_join(rsid, by = "id")
dim(myeloprof)
myeloprof <- myeloprof %>%
  select(id, rsids, alt, ref, af_alt, beta, sebeta, pval) %>%
  rename(ea_myeloprof = alt,
         oa_myeloprof = ref,
         eaf_myeloprof = af_alt,
         beta_myeloprof = beta,
         se_myeloprof = sebeta,
         p_myeloprof = pval) %>%
  select(id, rsids, ea_myeloprof, oa_myeloprof, eaf_myeloprof, beta_myeloprof, se_myeloprof, p_myeloprof)

myeloid <- myeloid %>%
  mutate(id = paste0(rsids, "_", pmin(ref, alt), "_", pmax(ref, alt))) %>%
  select(id, everything()) %>%
  inner_join(rsid, by = "id")
dim(myeloid)
myeloid <- myeloid %>%
  select(id, rsids, alt, ref, af_alt, beta, sebeta, pval) %>%
  rename(ea_myeloid = alt,
         oa_myeloid = ref,
         eaf_myeloid = af_alt,
         beta_myeloid = beta,
         se_myeloid = sebeta,
         p_myeloid = pval) %>%
  select(id, rsids, ea_myeloid, oa_myeloid, eaf_myeloid, beta_myeloid, se_myeloid, p_myeloid)

### Merge sffdr_novel2 with blood diseases:
ssfdr_blood <- sffdr_novel2 %>%
  left_join(aml, by = "id") %>%
  left_join(myelodys, by = "id") %>%
  left_join(myeloprof, by = "id") %>%
  left_join(myeloid, by = "id")
dim(ssfdr_blood)

### Harmonise data to CH risk increasing allele:
ea_mismatches <- ssfdr_blood %>%
  filter(ea != ea_aml |
           ea != ea_myelodys |
           ea != ea_myeloprof |
           ea != ea_myeloid)
ssfdr_blood_harmonised <- ssfdr_blood %>%
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
    
    beta_aml = ifelse(flip, -beta_aml, beta_aml),
    eaf_aml = ifelse(flip, 1 - eaf_aml, eaf_aml),
    or_aml = exp(beta_aml),
    l95cl_aml = exp(beta_aml - 1.96 * se_aml),
    u95cl_aml = exp(beta_aml + 1.96 * se_aml),
    
    beta_myelodys = ifelse(flip, -beta_myelodys, beta_myelodys),
    eaf_myelodys = ifelse(flip, 1 - eaf_myelodys, eaf_myelodys),
    or_myelodys = exp(beta_myelodys),
    l95cl_myelodys = exp(beta_myelodys - 1.96 * se_myelodys),
    u95cl_myelodys = exp(beta_myelodys + 1.96 * se_myelodys),
    
    beta_myeloprof = ifelse(flip, -beta_myeloprof, beta_myeloprof),
    eaf_myeloprof = ifelse(flip, 1 - eaf_myeloprof, eaf_myeloprof),
    or_myeloprof = exp(beta_myeloprof),
    l95cl_myeloprof = exp(beta_myeloprof - 1.96 * se_myeloprof),
    u95cl_myeloprof = exp(beta_myeloprof + 1.96 * se_myeloprof),
    
    beta_myeloid = ifelse(flip, -beta_myeloid, beta_myeloid),
    eaf_myeloid = ifelse(flip, 1 - eaf_myeloid, eaf_myeloid),
    or_myeloid = exp(beta_myeloid),
    l95cl_myeloid = exp(beta_myeloid - 1.96 * se_myeloid),
    u95cl_myeloid = exp(beta_myeloid + 1.96 * se_myeloid)) %>%
  select(-flip, -ea_orig, -oa_orig)

ssfdr_blood_harmonised  <- ssfdr_blood_harmonised  %>%
  select(-c(rsids.x, rsids.y, rsids.x.x, rsids.y.y,
            ea_aml, ea_myelodys, ea_myeloprof, ea_myeloid,
            oa_aml, oa_myelodys, oa_myeloprof, oa_myeloid)) %>%
  select(id, CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp,
         eaf_aml, beta_aml, se_aml, or_aml, l95cl_aml, u95cl_aml, p_aml,
         eaf_myelodys, beta_myelodys, se_myelodys, or_myelodys, l95cl_myelodys, u95cl_myelodys, p_myelodys,
         eaf_myeloprof, beta_myeloprof, se_myeloprof, or_myeloprof, l95cl_myeloprof, u95cl_myeloprof, p_myeloprof,
         eaf_myeloid, beta_myeloid, se_myeloid, or_myeloid, l95cl_myeloid, u95cl_myeloid, p_myeloid)

### Save:
write_csv(ssfdr_blood_harmonised, file = "sffdr_novel_finngen_aml_myelodys_myeloprof_myeloid.csv")