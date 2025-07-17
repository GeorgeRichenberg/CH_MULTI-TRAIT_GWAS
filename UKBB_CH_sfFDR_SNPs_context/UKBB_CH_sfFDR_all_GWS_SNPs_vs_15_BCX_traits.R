### Set wd and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
wbc_sffdr <- vroom("wbc_sffdr_harmonised.csv", col_select = -c(or, l95cl, u95cl, or, l95cl, u95cl)) 
neu_sffdr <- vroom("neu_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_neu", "oa_neu", "eaf_neu", "beta_neu", "se_neu", "pval_neu"))
lym_sffdr <- vroom("lym_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_lym", "oa_lym", "eaf_lym", "beta_lym", "se_lym", "pval_lym"))
mon_sffdr <- vroom("mon_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_mon", "oa_mon", "eaf_mon", "beta_mon", "se_mon", "pval_mon"))
bas_sffdr <- vroom("bas_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_bas", "oa_bas", "eaf_bas", "beta_bas", "se_bas", "pval_bas"))
eos_sffdr <- vroom("eos_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_eos", "oa_eos", "eaf_eos", "beta_eos", "se_eos", "pval_eos"))
rbc_sffdr <- vroom("rbc_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_rbc", "oa_rbc", "eaf_rbc", "beta_rbc", "se_rbc", "pval_rbc"))
hgb_sffdr <- vroom("hgb_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_hgb", "oa_hgb", "eaf_hgb", "beta_hgb", "se_hgb", "pval_hgb"))
hct_sffdr <- vroom("hct_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_hct", "oa_hct", "eaf_hct", "beta_hct", "se_hct", "pval_hct"))
mch_sffdr <- vroom("mch_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_mch", "oa_mch", "eaf_mch", "beta_mch", "se_mch", "pval_mch"))
mcv_sffdr <- vroom("mcv_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_mcv", "oa_mcv", "eaf_mcv", "beta_mcv", "se_mcv", "pval_mcv"))
mchc_sffdr <- vroom("mchc_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_mchc", "oa_mchc", "eaf_mchc", "beta_mchc", "se_mchc", "pval_mchc"))
rdw_sffdr <- vroom("rdw_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_rdw", "oa_rdw", "eaf_rdw", "beta_rdw", "se_rdw", "pval_rdw"))
plt_sffdr <- vroom("plt_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_plt", "oa_plt", "eaf_plt", "beta_plt", "se_plt", "pval_plt"))
mpv_sffdr <- vroom("mpv_sffdr_harmonised.csv", col_select = c("CH_GWAS", "rsID", "ea_mpv", "oa_mpv", "eaf_mpv", "beta_mpv", "se_mpv", "pval_mpv"))

### Create list:
datasets <- list(neu_sffdr,
                 lym_sffdr,
                 mon_sffdr,
                 bas_sffdr,
                 eos_sffdr,
                 rbc_sffdr,
                 hgb_sffdr,
                 hct_sffdr,
                 mch_sffdr,
                 mcv_sffdr,
                 mchc_sffdr,
                 rdw_sffdr,
                 plt_sffdr,
                 mpv_sffdr)

### Combine:
sffdr_joined <- wbc_sffdr
for (dataset in datasets) {
  sffdr_joined <- left_join(sffdr_joined, dataset, by = c("CH_GWAS", "rsID"))}
dim(sffdr_joined)

### Confirm each bcx trait is harmonised to sffdr:
ea_columns <- c("ea_wbc", "ea_neu", "ea_lym", "ea_mon", "ea_bas", "ea_eos", 
                "ea_rbc", "ea_hgb", "ea_hct", "ea_mch", "ea_mcv", "ea_mchc", "ea_rdw",
                "ea_plt", "ea_mpv")
diff_counts <- list()
for (col in ea_columns) {
  diff_counts[[col]] <- sum(sffdr_joined$ea != sffdr_joined[[col]])}
diff_counts

### Remove ea and oa columns:
sffdr_joined2 <- sffdr_joined %>%
  select(-c("ea_wbc", "ea_neu", "ea_lym", "ea_mon", "ea_bas", "ea_eos", 
            "ea_rbc", "ea_hgb", "ea_hct", "ea_mch", "ea_mcv", "ea_mchc", "ea_rdw", 
            "ea_plt", "ea_mpv",
            "oa_wbc", "oa_neu", "oa_lym", "oa_mon", "oa_bas", "oa_eos", 
            "oa_rbc", "oa_hgb", "oa_hct", "oa_mch", "oa_mcv", "oa_mchc", "oa_rdw", 
            "oa_plt", "oa_mpv"))
dim(sffdr_joined2)

### Harmonise to the risk increasing allele:
sffdr_joined3 <- sffdr_joined2 %>%
  mutate(flip = beta < 0,
         ea_orig = ea,
         oa_orig = oa,
         ea = ifelse(flip, oa_orig, ea),
         oa = ifelse(flip, ea_orig, oa),
         beta = ifelse(flip, -beta, beta),
         eaf = ifelse(flip, 1 - eaf, eaf),
         or = exp(beta),
         l95cl = exp(beta - 1.96 * se),
         u95cl = exp(beta + 1.96 * se)) %>%
  mutate(across(.cols = starts_with("beta_"),
                .fns = ~ ifelse(flip, -.x, .x),
                .names = "{.col}")) %>%
  mutate(across(.cols = starts_with("eaf_"),
                .fns = ~ ifelse(flip, 1 - .x, .x),
                .names = "{.col}")) %>%
  mutate(across(.cols = starts_with("beta_"),
                .fns = list(or = ~ exp(.x),
                            l95cl = ~ exp(.x - 1.96 * get(gsub("beta", "se", cur_column()))),
                            u95cl = ~ exp(.x + 1.96 * get(gsub("beta", "se", cur_column())))),
                            .names = "{.col}_{.fn}")) %>%
  select(-flip, -ea_orig, -oa_orig)

### Save:
write_csv(sffdr_joined3, file = "sffdr_all_gws_vs_bcx_traits_increasing_risk_allele.csv")