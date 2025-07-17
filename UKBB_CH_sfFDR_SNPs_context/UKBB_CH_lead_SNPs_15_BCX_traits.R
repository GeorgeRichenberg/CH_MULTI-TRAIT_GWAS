### Set wd and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)
library(MungeSumstats)

#### Load data:
kessler_lead <- vroom("kessler_lead_variants.csv", col_select = -c(Effect_UKB, LCI_Effect_UKB, UCI_Effect_UKB))
wbc <- vroom("WBC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_wbc"))
neu <- vroom("NEU_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_neu"))
lym <- vroom("LYM_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_lym"))
mon <- vroom("MON_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_mon"))
bas <- vroom("BAS_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_bas"))
eos <- vroom("EOS_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_eos"))
rbc <- vroom("RBC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_rbc"))
hgb <- vroom("HGB_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_hgb"))
hct <- vroom("HCT_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_hct"))
mch <- vroom("MCH_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_mch"))
mcv <- vroom("MCV_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_mcv"))
mchc <- vroom("MCHC_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_mchc"))
rdw <- vroom("RDW_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_rdw"))
plt <- vroom("PLT_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_plt"))
mpv <- vroom("MPV_Europeans_GWAMA_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz", col_select = -c(chr, bp, n_studies_bcx, n_samples_bcx)) %>%
  rename_with(~ str_replace(.x, "_bcx", "_mpv"))

### Filter lead variants from kessler_lead:
kessler_lead <- kessler_lead %>%
  group_by(locus) %>%
  slice_min(order_by = Pval_UKB, n = 1, with_ties = FALSE) %>%
  ungroup()

### Liftover from GRCh38 to GRCh37:
kessler_lead <- unite(kessler_lead, "SNP", Chr:Pos, remove = FALSE)
kessler_lead <- kessler_lead %>%
  rename(CHR = Chr,
         BP = Pos)
kessler_lead <- data.table(kessler_lead)
kessler_lead_lifted <- liftover(kessler_lead, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
kessler_lead_lifted <- as_tibble(kessler_lead_lifted)
dim(kessler_lead_lifted)
kessler_lead_lifted <- select(kessler_lead_lifted, -SNP)

### Create id column for kessler_lead variants:
kessler_lead_lifted <- kessler_lead_lifted %>%
  mutate(Allele1 = pmin(EA, OA),
         Allele2 = pmax(EA, OA),
         rs_number = paste0(CHR, ":", BP, "_", Allele1, "_", Allele2)) %>%
  select(rs_number, locus, everything()) %>%
  select(-c(ID, IMPUTATION_gen_build, Allele1, Allele2))

### Merge with each bcx trait:
kessler_lead_bcx <- kessler_lead_lifted %>%
  inner_join(wbc, by = "rs_number") %>%
  inner_join(neu, by = "rs_number") %>%
  inner_join(lym, by = "rs_number") %>%
  inner_join(mon, by = "rs_number") %>%
  inner_join(bas, by = "rs_number") %>%
  inner_join(eos, by = "rs_number") %>%
  inner_join(rbc, by = "rs_number") %>%
  inner_join(hgb, by = "rs_number") %>%
  inner_join(hct, by = "rs_number") %>%
  inner_join(mch, by = "rs_number") %>%
  inner_join(mcv, by = "rs_number") %>%
  inner_join(mchc, by = "rs_number") %>%
  inner_join(rdw, by = "rs_number") %>%
  inner_join(plt, by = "rs_number") %>%
  inner_join(mpv, by = "rs_number")

### Format with respect to the risk increasing allele:
kessler_harmonised <- kessler_lead_bcx %>%
  mutate(
    flip = Beta_UKB < 0,
    new_EA = if_else(flip, OA, EA),
    new_OA = if_else(flip, EA, OA),
    EAF_UKB = if_else(flip, 1 - EAF_UKB, EAF_UKB),
    Beta_UKB = abs(Beta_UKB)) %>%
  select(-c(EA, OA)) %>%
  rename(EA = new_EA, OA = new_OA) %>%
  select(rs_number, locus, CHR, BP, EA, OA, everything())

### Harmonise data:
harmonise_data <- function(data) {
  data %>%
    rowwise() %>%
    mutate(
      eaf_wbc = if_else(ea_wbc != EA, 1 - eaf_wbc, eaf_wbc),
      beta_wbc = if_else(ea_wbc != EA, -beta_wbc, beta_wbc),
      
      eaf_neu = if_else(ea_neu != EA, 1 - eaf_neu, eaf_neu),
      beta_neu = if_else(ea_neu != EA, -beta_neu, beta_neu),
      
      eaf_lym = if_else(ea_lym != EA, 1 - eaf_lym, eaf_lym),
      beta_lym = if_else(ea_lym != EA, -beta_lym, beta_lym),
      
      eaf_mon = if_else(ea_mon != EA, 1 - eaf_mon, eaf_mon),
      beta_mon = if_else(ea_mon != EA, -beta_mon, beta_mon),
      
      eaf_bas = if_else(ea_bas != EA, 1 - eaf_bas, eaf_bas),
      beta_bas = if_else(ea_bas != EA, -beta_bas, beta_bas),
      
      eaf_eos = if_else(ea_eos != EA, 1 - eaf_eos, eaf_eos),
      beta_eos = if_else(ea_eos != EA, -beta_eos, beta_eos),
      
      eaf_rbc = if_else(ea_rbc != EA, 1 - eaf_rbc, eaf_rbc),
      beta_rbc = if_else(ea_rbc != EA, -beta_rbc, beta_rbc),
      
      eaf_hgb = if_else(ea_hgb != EA, 1 - eaf_hgb, eaf_hgb),
      beta_hgb = if_else(ea_hgb != EA, -beta_hgb, beta_hgb),
      
      eaf_hct = if_else(ea_hct != EA, 1 - eaf_hct, eaf_hct),
      beta_hct = if_else(ea_hct != EA, -beta_hct, beta_hct),
      
      eaf_mch = if_else(ea_mch != EA, 1 - eaf_mch, eaf_mch),
      beta_mch = if_else(ea_mch != EA, -beta_mch, beta_mch),
      
      eaf_mcv = if_else(ea_mcv != EA, 1 - eaf_mcv, eaf_mcv),
      beta_mcv = if_else(ea_mcv != EA, -beta_mcv, beta_mcv),
      
      eaf_mchc = if_else(ea_mchc != EA, 1 - eaf_mchc, eaf_mchc),
      beta_mchc = if_else(ea_mchc != EA, -beta_mchc, beta_mchc),
      
      eaf_rdw = if_else(ea_rdw != EA, 1 - eaf_rdw, eaf_rdw),
      beta_rdw = if_else(ea_rdw != EA, -beta_rdw, beta_rdw),
      
      eaf_plt = if_else(ea_plt != EA, 1 - eaf_plt, eaf_plt),
      beta_plt = if_else(ea_plt != EA, -beta_plt, beta_plt),
      
      eaf_mpv = if_else(ea_mpv != EA, 1 - eaf_mpv, eaf_mpv),
      beta_mpv = if_else(ea_mpv != EA, -beta_mpv, beta_mpv)) %>%
    ungroup()}
kessler_harmonised2 <- harmonise_data(kessler_harmonised)
kessler_harmonised2
   
kessler_harmonised2 <- kessler_harmonised2 %>%
  select(-c(ea_wbc, ea_neu, ea_lym, ea_mon, ea_bas, ea_eos,
            ea_rbc, ea_hgb, ea_hct, ea_mch, ea_mcv, ea_mchc, ea_rdw, 
            ea_mpv, ea_plt,
            oa_wbc, oa_neu, oa_lym, oa_mon, oa_bas, oa_eos,
            oa_rbc, oa_hgb, oa_hct, oa_mch, oa_mcv, oa_mchc, oa_rdw, 
            oa_mpv, oa_plt))

### Calculate OR and CI:
kessler_harmonised3 <- kessler_harmonised2 %>%
  mutate(
    or_wbc = exp(beta_wbc),
    l95cl_wbc = exp(beta_wbc - 1.96 * se_wbc),
    u95cl_wbc = exp(beta_wbc + 1.96 * se_wbc),
    
    or_neu = exp(beta_neu),
    l95cl_neu = exp(beta_neu - 1.96 * se_neu),
    u95cl_neu = exp(beta_neu + 1.96 * se_neu),
    
    or_lym = exp(beta_lym),
    l95cl_lym = exp(beta_lym - 1.96 * se_lym),
    u95cl_lym = exp(beta_lym + 1.96 * se_lym),
    
    or_mon = exp(beta_mon),
    l95cl_mon = exp(beta_mon - 1.96 * se_mon),
    u95cl_mon = exp(beta_mon + 1.96 * se_mon),
    
    or_bas = exp(beta_bas),
    l95cl_bas = exp(beta_bas - 1.96 * se_bas),
    u95cl_bas = exp(beta_bas + 1.96 * se_bas),
    
    or_eos = exp(beta_eos),
    l95cl_eos = exp(beta_eos - 1.96 * se_eos),
    u95cl_eos = exp(beta_eos + 1.96 * se_eos),
    
    or_rbc = exp(beta_rbc),
    l95cl_rbc = exp(beta_rbc - 1.96 * se_rbc),
    u95cl_rbc = exp(beta_rbc + 1.96 * se_rbc),
    
    or_hgb = exp(beta_hgb),
    l95cl_hgb = exp(beta_hgb - 1.96 * se_hgb),
    u95cl_hgb = exp(beta_hgb + 1.96 * se_hgb),
    
    or_hct = exp(beta_hct),
    l95cl_hct = exp(beta_hct - 1.96 * se_hct),
    u95cl_hct = exp(beta_hct + 1.96 * se_hct),
    
    or_mch = exp(beta_mch),
    l95cl_mch = exp(beta_mch - 1.96 * se_mch),
    u95cl_mch = exp(beta_mch + 1.96 * se_mch),
    
    or_mcv = exp(beta_mcv),
    l95cl_mcv = exp(beta_mcv - 1.96 * se_mcv),
    u95cl_mcv = exp(beta_mcv + 1.96 * se_mcv),
    
    or_mchc = exp(beta_mchc),
    l95cl_mchc = exp(beta_mchc - 1.96 * se_mchc),
    u95cl_mchc = exp(beta_mchc + 1.96 * se_mchc),
    
    or_rdw = exp(beta_rdw),
    l95cl_rdw = exp(beta_rdw - 1.96 * se_rdw),
    u95cl_rdw = exp(beta_rdw + 1.96 * se_rdw),
    
    or_plt = exp(beta_plt),
    l95cl_plt = exp(beta_plt - 1.96 * se_plt),
    u95cl_plt = exp(beta_plt + 1.96 * se_plt),
    
    or_mpv = exp(beta_mpv),
    l95cl_mpv = exp(beta_mpv - 1.96 * se_mpv),
    u95cl_mpv = exp(beta_mpv + 1.96 * se_mpv))

### Order columns:
bcx_traits <- c("wbc", "neu", "lym", "mon", "bas", "eos",
                "rbc", "hgb", "hct", "mch", "mcv", "mchc", 
                "rdw", "plt", "mpv")
ordered_cols <- c(
  "rs_number", "locus", "CHR", "BP", "EA", "OA", "EAF_UKB", "Beta_UKB", "SE_UKB", "Pval_UKB", "nearestGene",
  unlist(lapply(bcx_traits, function(trait) {
    c(paste0("eaf_", trait),
      paste0("beta_", trait),
      paste0("se_", trait),
      paste0("or_", trait),
      paste0("l95cl_", trait),
      paste0("u95cl_", trait),
      paste0("pval_", trait))})))
kessler_harmonised3 <- kessler_harmonised3 %>% 
  select(all_of(ordered_cols))
kessler_harmonised3 <- kessler_harmonised3 %>%
  mutate(or_UKB = exp(Beta_UKB),
         l95cl_UKB = exp(Beta_UKB - 1.96 * SE_UKB),
         u95cl_UKB = exp(Beta_UKB + 1.96 * SE_UKB)) %>%
         select(rs_number, locus, CHR, BP, EA, OA, EAF_UKB, Beta_UKB, SE_UKB, or_UKB, l95cl_UKB, u95cl_UKB, Pval_UKB, nearestGene, everything())

### Provide rsID:
g1000 <- vroom("g1000_eur.bim", col_names = FALSE)
g1000 <- dplyr::select(g1000, -c(X3, X5, X6))
g1000 <- rename(g1000, rsid = X2)
g1000 <- unite(g1000, "SNP", c(X1, X4), remove = TRUE)
kessler_harmonised4 <- kessler_harmonised3 %>%
  mutate(SNP = paste(CHR, BP, sep = "_"))
  
kessler_harmonised4 <- kessler_harmonised4 %>%
  inner_join(g1000, by = "SNP")
kessler_harmonised4 <- kessler_harmonised4 %>%
  select(-rs_number) %>%
  select(rsid, everything())

### Save:
write_csv(kessler_harmonised4, file = "kessler_lead_snps_15_bcx_traits.csv")
  