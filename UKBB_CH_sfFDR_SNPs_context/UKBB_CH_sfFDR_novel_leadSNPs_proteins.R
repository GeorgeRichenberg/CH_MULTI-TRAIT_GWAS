### Set wd and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, p, fp))
cebpa <- vroom("cebpa_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))
flt3 <- vroom("flt3_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))
flt3lg <- vroom("flt3lg_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))
kit <- vroom("kit_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))
kitlg <- vroom("kitlg_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))
npm1 <- vroom("npm1_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))
tp53 <- vroom("tp53_all.csv", col_select = -c(SNP, CHR, BP, IMPUTATION_gen_build))

### Format gene data:
cebpa <- cebpa %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_cebpa = ALLELE1,
         oa_cebpa = ALLELE0,
         eaf_cebpa = A1FREQ,
         beta_cebpa = BETA,
         se_cebpa = SE,
         p_cebpa = p)
flt3 <- flt3 %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_flt3 = ALLELE1,
         oa_flt3 = ALLELE0,
         eaf_flt3 = A1FREQ,
         beta_flt3 = BETA,
         se_flt3 = SE,
         p_flt3 = p)
flt3lg <- flt3lg %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_flt3lg = ALLELE1,
         oa_flt3lg = ALLELE0,
         eaf_flt3lg = A1FREQ,
         beta_flt3lg = BETA,
         se_flt3lg = SE,
         p_flt3lg = p)
kit <- kit %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_kit = ALLELE1,
         oa_kit = ALLELE0,
         eaf_kit = A1FREQ,
         beta_kit = BETA,
         se_kit = SE,
         p_kit = p)
kitlg <- kitlg %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_kitlg = ALLELE1,
         oa_kitlg = ALLELE0,
         eaf_kitlg = A1FREQ,
         beta_kitlg = BETA,
         se_kitlg = SE,
         p_kitlg = p)
npm1 <- npm1 %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_npm1 = ALLELE1,
         oa_npm1 = ALLELE0, 
         eaf_npm1 = A1FREQ,
         beta_npm1 = BETA,
         se_npm1 = SE,
         p_npm1 = p)
tp53 <- tp53 %>%
  mutate(p = 10^-LOG10P) %>%
  select(-LOG10P) %>%
  rename(rsID = rsid,
         ea_tp53 = ALLELE1,
         oa_tp53 = ALLELE0, 
         eaf_tp53 = A1FREQ,
         beta_tp53 = BETA,
         se_tp53 = SE,
         p_tp53 = p)

### Merge data:
sffdr_novel <- sffdr_novel %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR|\\[|\\]"))
sffdr_gene <- sffdr_novel %>%
  inner_join(cebpa, by = "rsID") %>%
  inner_join(flt3, by = "rsID") %>%
  inner_join(flt3lg, by = "rsID") %>%
  inner_join(kit, by = "rsID") %>%
  inner_join(kitlg, by = "rsID") %>%
  inner_join(npm1, by = "rsID") %>%
  inner_join(tp53, by = "rsID")

### Harmonise data to CH risk increasing allele:
ea_mismatches <- sffdr_gene %>%
  filter(ea != ea_cebpa |
         ea != ea_flt3 |
         ea != ea_flt3lg |
         ea != ea_kit,
         ea != ea_kitlg,
         ea != ea_npm1,
         ea != ea_tp53)
ea_mismatches
sffdr_gene <- sffdr_gene %>%
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
    
    beta_cebpa = ifelse(flip, -beta_cebpa, beta_cebpa),
    eaf_cebpa = ifelse(flip, 1 - eaf_cebpa, eaf_cebpa),
    or_cebpa = exp(beta_cebpa),
    l95cl_cebpa = exp(beta_cebpa - 1.96 * se_cebpa),
    u95cl_cebpa = exp(beta_cebpa + 1.96 * se_cebpa),
    
    beta_flt3 = ifelse(flip, -beta_flt3, beta_flt3),
    eaf_flt3 = ifelse(flip, 1 - eaf_flt3, eaf_flt3),
    or_flt3 = exp(beta_flt3),
    l95cl_flt3 = exp(beta_flt3 - 1.96 * se_flt3),
    u95cl_flt3 = exp(beta_flt3 + 1.96 * se_flt3),
    
    beta_flt3lg = ifelse(flip, -beta_flt3lg, beta_flt3lg),
    eaf_flt3lg = ifelse(flip, 1 - eaf_flt3lg, eaf_flt3lg),
    or_flt3lg = exp(beta_flt3lg),
    l95cl_flt3lg = exp(beta_flt3lg - 1.96 * se_flt3lg),
    u95cl_flt3lg = exp(beta_flt3lg + 1.96 * se_flt3lg),
    
    beta_kit = ifelse(flip, -beta_kit, beta_kit),
    eaf_kit = ifelse(flip, 1 - eaf_kit, eaf_kit),
    or_kit = exp(beta_kit),
    l95cl_kit = exp(beta_kit - 1.96 * se_kit),
    u95cl_kit = exp(beta_kit + 1.96 * se_kit),
    
    beta_kitlg = ifelse(flip, -beta_kitlg, beta_kitlg),
    eaf_kitlg = ifelse(flip, 1 - eaf_kitlg, eaf_kitlg),
    or_kitlg = exp(beta_kitlg),
    l95cl_kitlg = exp(beta_kitlg - 1.96 * se_kitlg),
    u95cl_kitlg = exp(beta_kitlg + 1.96 * se_kitlg),
    
    beta_npm1 = ifelse(flip, -beta_npm1, beta_npm1),
    eaf_npm1 = ifelse(flip, 1 - eaf_npm1, eaf_npm1),
    or_npm1 = exp(beta_npm1),
    l95cl_npm1 = exp(beta_npm1 - 1.96 * se_npm1),
    u95cl_npm1 = exp(beta_npm1 + 1.96 * se_npm1),
    
    beta_tp53 = ifelse(flip, -beta_tp53, beta_tp53),
    eaf_tp53 = ifelse(flip, 1 - eaf_tp53, eaf_tp53),
    or_tp53 = exp(beta_tp53),
    l95cl_tp53 = exp(beta_tp53 - 1.96 * se_tp53),
    u95cl_tp53 = exp(beta_tp53 + 1.96 * se_tp53)) %>%
  select(-flip, -ea_orig, -oa_orig)
sffdr_gene

### Format and save:
sffdr_gene <- sffdr_gene %>%
  select(-c(ea_cebpa, ea_flt3, ea_flt3lg, ea_kit, ea_kitlg, ea_npm1, ea_tp53,
            oa_cebpa, oa_flt3, oa_flt3lg, oa_kit, oa_kitlg, oa_npm1, oa_tp53)) %>%
  select(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, 
         eaf_cebpa, beta_cebpa, se_cebpa, or_cebpa, l95cl_cebpa, u95cl_cebpa, p_cebpa, 
         eaf_flt3, beta_flt3, se_flt3, or_flt3, l95cl_flt3, u95cl_flt3, p_flt3, 
         eaf_flt3lg, beta_flt3lg, se_flt3lg, or_flt3lg, l95cl_flt3lg, u95cl_flt3lg, p_flt3lg, 
         eaf_kit, beta_kit, se_kit, or_kit, l95cl_kit, u95cl_kit, p_kit, 
         eaf_kitlg, beta_kitlg, se_kitlg, or_kitlg, l95cl_kitlg, u95cl_kitlg, p_kitlg, 
         eaf_npm1, beta_npm1, se_npm1, or_npm1, l95cl_npm1, u95cl_npm1, p_npm1, 
         eaf_tp53, beta_tp53, se_tp53, or_tp53, l95cl_tp53, u95cl_tp53, p_tp53)
write_csv(sffdr_gene, file = "sffdr_novel_plasma_proteins_cebpa_flt3_flt3lg_kit_kitlg_npm1_tp53.csv")