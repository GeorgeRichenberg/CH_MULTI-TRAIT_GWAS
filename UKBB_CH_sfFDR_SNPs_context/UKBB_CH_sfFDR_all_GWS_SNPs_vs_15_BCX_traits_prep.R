### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
sffdr <- vroom("sffdr_all_ch_gws_snps_st1.csv")

### Format sffdr:
sffdr2 <- sffdr %>%
  mutate(alleles = pmap_chr(list(ea, oa), ~paste0(sort(c(..1, ..2)), collapse = "_")),
         rs_number = paste0(chr, ":", pos, "_", alleles)) %>%
  select(-alleles)
rs_number <- sffdr2 %>%
  select(CH_GWAS, rs_number)

### List BCX types:
bcx_types <- c("WBC", "NEU", "LYM", "MON", "BAS", "EOS", 
               "RBC", "HGB", "HCT", "MCH", "MCV", "MCHC", "RDW", 
               "PLT", "MPV")

### Create function to format bcx data:
format_bcx <- function(bcx) {
  df <- vroom(paste0(bcx, "_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz"),
              col_select = c(rs_number, ea_bcx, oa_bcx, eaf_bcx, beta_bcx, se_bcx, pval_bcx))
  renamed_cols <- c("ea", "oa", "eaf", "beta", "se", "pval")
  names(df)[2:7] <- paste0(renamed_cols, "_", tolower(bcx))
  df <- inner_join(rs_number, df, by = "rs_number")
  return(df)}

### Format bcx data:
wbc_sffdr  <- format_bcx("WBC")
dim(wbc_sffdr)
rm(wbc)

neu_sffdr  <- format_bcx("NEU")
dim(neu_sffdr)
rm(neu)

lym_sffdr  <- format_bcx("LYM")
dim(lym_sffdr)
rm(lym)

mon_sffdr  <- format_bcx("MON")
dim(mon_sffdr)
rm(mon)

bas_sffdr  <- format_bcx("BAS")
dim(bas_sffdr)
rm(bas)

eos_sffdr  <- format_bcx("EOS")
dim(eos_sffdr)
rm(eos)

rbc_sffdr  <- format_bcx("RBC")
dim(rbc_sffdr)
rm(rbc)

hgb_sffdr  <- format_bcx("HGB")
dim(hgb_sffdr)
rm(hgb)

hct_sffdr  <- format_bcx("HCT")
dim(hct_sffdr)
rm(hct)

mch_sffdr  <- format_bcx("MCH")
dim(mch_sffdr)
rm(mch)

mcv_sffdr  <- format_bcx("MCV")
dim(mcv_sffdr)
rm(mcv)

mchc_sffdr <- format_bcx("MCHC")
dim(mchc_sffdr)
rm(mchc)

rdw_sffdr  <- format_bcx("RDW")
dim(rdw_sffdr)
rm(rdw)

plt_sffdr  <- format_bcx("PLT")
dim(plt_sffdr)
rm(plt)

mpv_sffdr  <- format_bcx("MPV")
dim(mpv_sffdr)
rm(mpv)

### Harmonise each dataset with sffdr:
wbc_sffdr <- sffdr2 %>%
  inner_join(wbc_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(ea_wbc_harmonised = case_when(ea == ea_wbc ~ ea_wbc,  
                                       ea == oa_wbc ~ oa_wbc),
         oa_wbc_harmonised = case_when(ea == ea_wbc ~ oa_wbc,  
                                       ea == oa_wbc ~ ea_wbc),
         eaf_wbc_harmonised = case_when(ea == ea_wbc ~ eaf_wbc,     
                                        ea == oa_wbc ~ 1 - eaf_wbc),
         beta_wbc_harmonised = case_when(ea == ea_wbc ~ beta_wbc,  
                                         ea == oa_wbc ~ -beta_wbc)) %>%
  select(-c(ea_wbc, oa_wbc, eaf_wbc, beta_wbc)) %>%
  mutate(or_wbc = exp(beta_wbc_harmonised),  
         l95cl_wbc = exp(beta_wbc_harmonised - (1.96 * se_wbc)),  
         u95cl_wbc = exp(beta_wbc_harmonised + (1.96 * se_wbc))) %>%
  rename(ea_wbc = ea_wbc_harmonised,
         oa_wbc = oa_wbc_harmonised,
         eaf_wbc = eaf_wbc_harmonised,
         beta_wbc = beta_wbc_harmonised) %>%
  select(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
         ea_wbc, oa_wbc, eaf_wbc, beta_wbc, se_wbc, or_wbc, l95cl_wbc, u95cl_wbc, pval_wbc)
dim(wbc_sffdr)

neu_sffdr <- sffdr2 %>%
  inner_join(neu_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(ea_neu_harmonised = case_when(ea == ea_neu ~ ea_neu,  
                                       ea == oa_neu ~ oa_neu),
         oa_neu_harmonised = case_when(ea == ea_neu ~ oa_neu,  
                                       ea == oa_neu ~ ea_neu),
         eaf_neu_harmonised = case_when(ea == ea_neu ~ eaf_neu,     
                                        ea == oa_neu ~ 1 - eaf_neu),
         beta_neu_harmonised = case_when(ea == ea_neu ~ beta_neu,  
                                         ea == oa_neu ~ -beta_neu)) %>%
  select(-c(ea_neu, oa_neu, eaf_neu, beta_neu)) %>%
  mutate(or_neu = exp(beta_neu_harmonised),  
         l95cl_neu = exp(beta_neu_harmonised - (1.96 * se_neu)),  
         u95cl_neu = exp(beta_neu_harmonised + (1.96 * se_neu))) %>%
  rename(ea_neu = ea_neu_harmonised,
         oa_neu = oa_neu_harmonised,
         eaf_neu = eaf_neu_harmonised,
         beta_neu = beta_neu_harmonised) %>%
  select(CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
         ea_neu, oa_neu, eaf_neu, beta_neu, se_neu, or_neu, l95cl_neu, u95cl_neu, pval_neu)
dim(neu_sffdr)

lym_sffdr <- sffdr2 %>%
  inner_join(lym_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_lym_harmonised = case_when(
      ea == ea_lym ~ ea_lym,  
      ea == oa_lym ~ oa_lym),
    oa_lym_harmonised = case_when(
      ea == ea_lym ~ oa_lym,  
      ea == oa_lym ~ ea_lym),
    eaf_lym_harmonised = case_when(
      ea == ea_lym ~ eaf_lym,     
      ea == oa_lym ~ 1 - eaf_lym),
    beta_lym_harmonised = case_when(
      ea == ea_lym ~ beta_lym,  
      ea == oa_lym ~ -beta_lym)) %>%
  select(-c(ea_lym, oa_lym, eaf_lym, beta_lym)) %>%
  mutate(
    or_lym = exp(beta_lym_harmonised),  
    l95cl_lym = exp(beta_lym_harmonised - (1.96 * se_lym)),  
    u95cl_lym = exp(beta_lym_harmonised + (1.96 * se_lym))) %>%
  rename(
    ea_lym = ea_lym_harmonised,
    oa_lym = oa_lym_harmonised,
    eaf_lym = eaf_lym_harmonised,
    beta_lym = beta_lym_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_lym, oa_lym, eaf_lym, beta_lym, se_lym, or_lym, l95cl_lym, u95cl_lym, pval_lym)
dim(lym_sffdr)

mon_sffdr <- sffdr2 %>%
  inner_join(mon_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_mon_harmonised = case_when(
      ea == ea_mon ~ ea_mon,  
      ea == oa_mon ~ oa_mon),
    oa_mon_harmonised = case_when(
      ea == ea_mon ~ oa_mon,  
      ea == oa_mon ~ ea_mon),
    eaf_mon_harmonised = case_when(
      ea == ea_mon ~ eaf_mon,     
      ea == oa_mon ~ 1 - eaf_mon),
    beta_mon_harmonised = case_when(
      ea == ea_mon ~ beta_mon,  
      ea == oa_mon ~ -beta_mon)) %>%
  select(-c(ea_mon, oa_mon, eaf_mon, beta_mon)) %>%
  mutate(
    or_mon = exp(beta_mon_harmonised),  
    l95cl_mon = exp(beta_mon_harmonised - (1.96 * se_mon)),  
    u95cl_mon = exp(beta_mon_harmonised + (1.96 * se_mon))) %>%
  rename(
    ea_mon = ea_mon_harmonised,
    oa_mon = oa_mon_harmonised,
    eaf_mon = eaf_mon_harmonised,
    beta_mon = beta_mon_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_mon, oa_mon, eaf_mon, beta_mon, se_mon, or_mon, l95cl_mon, u95cl_mon, pval_mon)
dim(mon_sffdr)

bas_sffdr <- sffdr2 %>%
  inner_join(bas_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_bas_harmonised = case_when(
      ea == ea_bas ~ ea_bas,  
      ea == oa_bas ~ oa_bas),
    oa_bas_harmonised = case_when(
      ea == ea_bas ~ oa_bas,  
      ea == oa_bas ~ ea_bas),
    eaf_bas_harmonised = case_when(
      ea == ea_bas ~ eaf_bas,     
      ea == oa_bas ~ 1 - eaf_bas),
    beta_bas_harmonised = case_when(
      ea == ea_bas ~ beta_bas,  
      ea == oa_bas ~ -beta_bas)) %>%
  select(-c(ea_bas, oa_bas, eaf_bas, beta_bas)) %>%
  mutate(
    or_bas = exp(beta_bas_harmonised),  
    l95cl_bas = exp(beta_bas_harmonised - (1.96 * se_bas)),  
    u95cl_bas = exp(beta_bas_harmonised + (1.96 * se_bas))) %>%
  rename(
    ea_bas = ea_bas_harmonised,
    oa_bas = oa_bas_harmonised,
    eaf_bas = eaf_bas_harmonised,
    beta_bas = beta_bas_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_bas, oa_bas, eaf_bas, beta_bas, se_bas, or_bas, l95cl_bas, u95cl_bas, pval_bas)
dim(bas_sffdr)

eos_sffdr <- sffdr2 %>%
  inner_join(eos_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_eos_harmonised = case_when(
      ea == ea_eos ~ ea_eos,  
      ea == oa_eos ~ oa_eos),
    oa_eos_harmonised = case_when(
      ea == ea_eos ~ oa_eos,  
      ea == oa_eos ~ ea_eos),
    eaf_eos_harmonised = case_when(
      ea == ea_eos ~ eaf_eos,     
      ea == oa_eos ~ 1 - eaf_eos),
    beta_eos_harmonised = case_when(
      ea == ea_eos ~ beta_eos,  
      ea == oa_eos ~ -beta_eos)) %>%
  select(-c(ea_eos, oa_eos, eaf_eos, beta_eos)) %>%
  mutate(
    or_eos = exp(beta_eos_harmonised),  
    l95cl_eos = exp(beta_eos_harmonised - (1.96 * se_eos)),  
    u95cl_eos = exp(beta_eos_harmonised + (1.96 * se_eos))) %>%
  rename(
    ea_eos = ea_eos_harmonised,
    oa_eos = oa_eos_harmonised,
    eaf_eos = eaf_eos_harmonised,
    beta_eos = beta_eos_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_eos, oa_eos, eaf_eos, beta_eos, se_eos, or_eos, l95cl_eos, u95cl_eos, pval_eos)
dim(eos_sffdr)

rbc_sffdr <- sffdr2 %>%
  inner_join(rbc_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_rbc_harmonised = case_when(
      ea == ea_rbc ~ ea_rbc,  
      ea == oa_rbc ~ oa_rbc),
    oa_rbc_harmonised = case_when(
      ea == ea_rbc ~ oa_rbc,  
      ea == oa_rbc ~ ea_rbc),
    eaf_rbc_harmonised = case_when(
      ea == ea_rbc ~ eaf_rbc,     
      ea == oa_rbc ~ 1 - eaf_rbc),
    beta_rbc_harmonised = case_when(
      ea == ea_rbc ~ beta_rbc,  
      ea == oa_rbc ~ -beta_rbc)) %>%
  select(-c(ea_rbc, oa_rbc, eaf_rbc, beta_rbc)) %>%
  mutate(
    or_rbc = exp(beta_rbc_harmonised),  
    l95cl_rbc = exp(beta_rbc_harmonised - (1.96 * se_rbc)),  
    u95cl_rbc = exp(beta_rbc_harmonised + (1.96 * se_rbc))) %>%
  rename(
    ea_rbc = ea_rbc_harmonised,
    oa_rbc = oa_rbc_harmonised,
    eaf_rbc = eaf_rbc_harmonised,
    beta_rbc = beta_rbc_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_rbc, oa_rbc, eaf_rbc, beta_rbc, se_rbc, or_rbc, l95cl_rbc, u95cl_rbc, pval_rbc)
dim(rbc_sffdr)

hgb_sffdr <- sffdr2 %>%
  inner_join(hgb_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_hgb_harmonised = case_when(
      ea == ea_hgb ~ ea_hgb,  
      ea == oa_hgb ~ oa_hgb),
    oa_hgb_harmonised = case_when(
      ea == ea_hgb ~ oa_hgb,  
      ea == oa_hgb ~ ea_hgb),
    eaf_hgb_harmonised = case_when(
      ea == ea_hgb ~ eaf_hgb,     
      ea == oa_hgb ~ 1 - eaf_hgb),
    beta_hgb_harmonised = case_when(
      ea == ea_hgb ~ beta_hgb,  
      ea == oa_hgb ~ -beta_hgb)) %>%
  select(-c(ea_hgb, oa_hgb, eaf_hgb, beta_hgb)) %>%
  mutate(
    or_hgb = exp(beta_hgb_harmonised),  
    l95cl_hgb = exp(beta_hgb_harmonised - (1.96 * se_hgb)),  
    u95cl_hgb = exp(beta_hgb_harmonised + (1.96 * se_hgb))) %>%
  rename(
    ea_hgb = ea_hgb_harmonised,
    oa_hgb = oa_hgb_harmonised,
    eaf_hgb = eaf_hgb_harmonised,
    beta_hgb = beta_hgb_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_hgb, oa_hgb, eaf_hgb, beta_hgb, se_hgb, or_hgb, l95cl_hgb, u95cl_hgb, pval_hgb)
dim(hgb_sffdr)

hct_sffdr <- sffdr2 %>%
  inner_join(hct_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_hct_harmonised = case_when(
      ea == ea_hct ~ ea_hct,  
      ea == oa_hct ~ oa_hct),
    oa_hct_harmonised = case_when(
      ea == ea_hct ~ oa_hct,  
      ea == oa_hct ~ ea_hct),
    eaf_hct_harmonised = case_when(
      ea == ea_hct ~ eaf_hct,     
      ea == oa_hct ~ 1 - eaf_hct),
    beta_hct_harmonised = case_when(
      ea == ea_hct ~ beta_hct,  
      ea == oa_hct ~ -beta_hct)) %>%
  select(-c(ea_hct, oa_hct, eaf_hct, beta_hct)) %>%
  mutate(
    or_hct = exp(beta_hct_harmonised),  
    l95cl_hct = exp(beta_hct_harmonised - (1.96 * se_hct)),  
    u95cl_hct = exp(beta_hct_harmonised + (1.96 * se_hct))) %>%
  rename(
    ea_hct = ea_hct_harmonised,
    oa_hct = oa_hct_harmonised,
    eaf_hct = eaf_hct_harmonised,
    beta_hct = beta_hct_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_hct, oa_hct, eaf_hct, beta_hct, se_hct, or_hct, l95cl_hct, u95cl_hct, pval_hct)
dim(hct_sffdr)

mch_sffdr <- sffdr2 %>%
  inner_join(mch_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_mch_harmonised = case_when(
      ea == ea_mch ~ ea_mch,  
      ea == oa_mch ~ oa_mch),
    oa_mch_harmonised = case_when(
      ea == ea_mch ~ oa_mch,  
      ea == oa_mch ~ ea_mch),
    eaf_mch_harmonised = case_when(
      ea == ea_mch ~ eaf_mch,     
      ea == oa_mch ~ 1 - eaf_mch),
    beta_mch_harmonised = case_when(
      ea == ea_mch ~ beta_mch,  
      ea == oa_mch ~ -beta_mch)) %>%
  select(-c(ea_mch, oa_mch, eaf_mch, beta_mch)) %>%
  mutate(
    or_mch = exp(beta_mch_harmonised),  
    l95cl_mch = exp(beta_mch_harmonised - (1.96 * se_mch)),  
    u95cl_mch = exp(beta_mch_harmonised + (1.96 * se_mch))) %>%
  rename(
    ea_mch = ea_mch_harmonised,
    oa_mch = oa_mch_harmonised,
    eaf_mch = eaf_mch_harmonised,
    beta_mch = beta_mch_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_mch, oa_mch, eaf_mch, beta_mch, se_mch, or_mch, l95cl_mch, u95cl_mch, pval_mch)
dim(mch_sffdr)

mcv_sffdr <- sffdr2 %>%
  inner_join(mcv_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_mcv_harmonised = case_when(
      ea == ea_mcv ~ ea_mcv,  
      ea == oa_mcv ~ oa_mcv),
    oa_mcv_harmonised = case_when(
      ea == ea_mcv ~ oa_mcv,  
      ea == oa_mcv ~ ea_mcv),
    eaf_mcv_harmonised = case_when(
      ea == ea_mcv ~ eaf_mcv,     
      ea == oa_mcv ~ 1 - eaf_mcv),
    beta_mcv_harmonised = case_when(
      ea == ea_mcv ~ beta_mcv,  
      ea == oa_mcv ~ -beta_mcv)) %>%
  select(-c(ea_mcv, oa_mcv, eaf_mcv, beta_mcv)) %>%
  mutate(
    or_mcv = exp(beta_mcv_harmonised),  
    l95cl_mcv = exp(beta_mcv_harmonised - (1.96 * se_mcv)),  
    u95cl_mcv = exp(beta_mcv_harmonised + (1.96 * se_mcv))) %>%
  rename(
    ea_mcv = ea_mcv_harmonised,
    oa_mcv = oa_mcv_harmonised,
    eaf_mcv = eaf_mcv_harmonised,
    beta_mcv = beta_mcv_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_mcv, oa_mcv, eaf_mcv, beta_mcv, se_mcv, or_mcv, l95cl_mcv, u95cl_mcv, pval_mcv)
dim(mcv_sffdr)

mchc_sffdr <- sffdr2 %>%
  inner_join(mchc_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_mchc_harmonised = case_when(
      ea == ea_mchc ~ ea_mchc,  
      ea == oa_mchc ~ oa_mchc),
    oa_mchc_harmonised = case_when(
      ea == ea_mchc ~ oa_mchc,  
      ea == oa_mchc ~ ea_mchc),
    eaf_mchc_harmonised = case_when(
      ea == ea_mchc ~ eaf_mchc,     
      ea == oa_mchc ~ 1 - eaf_mchc),
    beta_mchc_harmonised = case_when(
      ea == ea_mchc ~ beta_mchc,  
      ea == oa_mchc ~ -beta_mchc)) %>%
  select(-c(ea_mchc, oa_mchc, eaf_mchc, beta_mchc)) %>%
  mutate(
    or_mchc = exp(beta_mchc_harmonised),  
    l95cl_mchc = exp(beta_mchc_harmonised - (1.96 * se_mchc)),  
    u95cl_mchc = exp(beta_mchc_harmonised + (1.96 * se_mchc))) %>%
  rename(
    ea_mchc = ea_mchc_harmonised,
    oa_mchc = oa_mchc_harmonised,
    eaf_mchc = eaf_mchc_harmonised,
    beta_mchc = beta_mchc_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_mchc, oa_mchc, eaf_mchc, beta_mchc, se_mchc, or_mchc, l95cl_mchc, u95cl_mchc, pval_mchc)
dim(mchc_sffdr)

rdw_sffdr <- sffdr2 %>%
  inner_join(rdw_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_rdw_harmonised = case_when(
      ea == ea_rdw ~ ea_rdw,  
      ea == oa_rdw ~ oa_rdw),
    oa_rdw_harmonised = case_when(
      ea == ea_rdw ~ oa_rdw,  
      ea == oa_rdw ~ ea_rdw),
    eaf_rdw_harmonised = case_when(
      ea == ea_rdw ~ eaf_rdw,     
      ea == oa_rdw ~ 1 - eaf_rdw),
    beta_rdw_harmonised = case_when(
      ea == ea_rdw ~ beta_rdw,  
      ea == oa_rdw ~ -beta_rdw)) %>%
  select(-c(ea_rdw, oa_rdw, eaf_rdw, beta_rdw)) %>%
  mutate(
    or_rdw = exp(beta_rdw_harmonised),  
    l95cl_rdw = exp(beta_rdw_harmonised - (1.96 * se_rdw)),  
    u95cl_rdw = exp(beta_rdw_harmonised + (1.96 * se_rdw))) %>%
  rename(
    ea_rdw = ea_rdw_harmonised,
    oa_rdw = oa_rdw_harmonised,
    eaf_rdw = eaf_rdw_harmonised,
    beta_rdw = beta_rdw_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_rdw, oa_rdw, eaf_rdw, beta_rdw, se_rdw, or_rdw, l95cl_rdw, u95cl_rdw, pval_rdw)
dim(rdw_sffdr)

plt_sffdr <- sffdr2 %>%
  inner_join(plt_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_plt_harmonised = case_when(
      ea == ea_plt ~ ea_plt,  
      ea == oa_plt ~ oa_plt),
    oa_plt_harmonised = case_when(
      ea == ea_plt ~ oa_plt,  
      ea == oa_plt ~ ea_plt),
    eaf_plt_harmonised = case_when(
      ea == ea_plt ~ eaf_plt,     
      ea == oa_plt ~ 1 - eaf_plt),
    beta_plt_harmonised = case_when(
      ea == ea_plt ~ beta_plt,  
      ea == oa_plt ~ -beta_plt)) %>%
  select(-c(ea_plt, oa_plt, eaf_plt, beta_plt)) %>%
  mutate(
    or_plt = exp(beta_plt_harmonised),  
    l95cl_plt = exp(beta_plt_harmonised - (1.96 * se_plt)),  
    u95cl_plt = exp(beta_plt_harmonised + (1.96 * se_plt))) %>%
  rename(
    ea_plt = ea_plt_harmonised,
    oa_plt = oa_plt_harmonised,
    eaf_plt = eaf_plt_harmonised,
    beta_plt = beta_plt_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_plt, oa_plt, eaf_plt, beta_plt, se_plt, or_plt, l95cl_plt, u95cl_plt, pval_plt)
dim(plt_sffdr)

mpv_sffdr <- sffdr2 %>%
  inner_join(mpv_sffdr, by = c("rs_number", "CH_GWAS")) %>%
  mutate(
    ea_mpv_harmonised = case_when(
      ea == ea_mpv ~ ea_mpv,  
      ea == oa_mpv ~ oa_mpv),
    oa_mpv_harmonised = case_when(
      ea == ea_mpv ~ oa_mpv,  
      ea == oa_mpv ~ ea_mpv),
    eaf_mpv_harmonised = case_when(
      ea == ea_mpv ~ eaf_mpv,     
      ea == oa_mpv ~ 1 - eaf_mpv),
    beta_mpv_harmonised = case_when(
      ea == ea_mpv ~ beta_mpv,  
      ea == oa_mpv ~ -beta_mpv)) %>%
  select(-c(ea_mpv, oa_mpv, eaf_mpv, beta_mpv)) %>%
  mutate(
    or_mpv = exp(beta_mpv_harmonised),  
    l95cl_mpv = exp(beta_mpv_harmonised - (1.96 * se_mpv)),  
    u95cl_mpv = exp(beta_mpv_harmonised + (1.96 * se_mpv))) %>%
  rename(
    ea_mpv = ea_mpv_harmonised,
    oa_mpv = oa_mpv_harmonised,
    eaf_mpv = eaf_mpv_harmonised,
    beta_mpv = beta_mpv_harmonised) %>%
  select(
    CH_GWAS, rsID, chr, pos, ea, oa, eaf, beta, se, or, l95cl, u95cl, p, fp, nearestGene, 
    ea_mpv, oa_mpv, eaf_mpv, beta_mpv, se_mpv, or_mpv, l95cl_mpv, u95cl_mpv, pval_mpv)
dim(mpv_sffdr)

### Save
write_csv(wbc_sffdr, "wbc_sffdr_harmonised.csv")
write_csv(neu_sffdr, "neu_sffdr_harmonised.csv")
write_csv(lym_sffdr, "lym_sffdr_harmonised.csv")
write_csv(mon_sffdr, "mon_sffdr_harmonised.csv")
write_csv(bas_sffdr, "bas_sffdr_harmonised.csv")
write_csv(eos_sffdr, "eos_sffdr_harmonised.csv")
write_csv(rbc_sffdr, "rbc_sffdr_harmonised.csv")
write_csv(hgb_sffdr, "hgb_sffdr_harmonised.csv")
write_csv(hct_sffdr, "hct_sffdr_harmonised.csv")
write_csv(mch_sffdr, "mch_sffdr_harmonised.csv")
write_csv(mcv_sffdr, "mcv_sffdr_harmonised.csv")
write_csv(mchc_sffdr, "mchc_sffdr_harmonised.csv")
write_csv(rdw_sffdr, "rdw_sffdr_harmonised.csv")
write_csv(plt_sffdr, "plt_sffdr_harmonised.csv")
write_csv(mpv_sffdr, "mpv_sffdr_harmonised.csv")