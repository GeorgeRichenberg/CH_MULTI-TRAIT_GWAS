### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)
library(MungeSumstats)

### Load sfFDR novel SNPs: 
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv")

### Lift sffdr_novel to GRCh38:
sffdr_novel <- sffdr_novel %>%
  rename(CHR = chr,
         BP = pos)
sffdr_novel <- unite(sffdr_novel, "SNP", CHR:BP, remove = FALSE)
sffdr_novel <- data.table(sffdr_novel)
sffdr_novel_lifted <- liftover(sffdr_novel, convert_ref_genome = "GRCh38", ref_genome = "GRCh37")
sffdr_novel_lifted <- as_tibble(sffdr_novel_lifted)
sffdr_novel_lifted <- select(sffdr_novel_lifted, -SNP)
sffdr_novel_lifted <- unite(sffdr_novel_lifted, "SNP", CHR:BP, remove = FALSE)

### Define target SNPs:
target_snps <- tribble(
  ~rsid,          ~chr, ~pos,
  "rs34016663",    1, 26859634,
  "rs1707334",     1, 46140056,
  "rs11674811",    2, 37011023,
  "rs113542380",   2, 43237679,
  "rs342467",      4, 87131067,
  "rs1062158",     5, 142143435,
  "rs2237154",     6, 15399883,
  "rs9404952",     6, 29836388,
  "rs536897754",   6, 30451389,
  "rs1054684",     6, 32809144,
  "rs35150201",    7, 135661514,
  "rs55903902",   11, 61764519,
  "rs2417347",    12, 14404839,
  "rs60626639",   16, 67591894,
  "rs4262954",    16, 87780989,
  "rs17758695",   18, 63253621,
  "rs1474547",    21, 14981960,
  "rs1034332",    21, 33976984,
  "rs9267820",     6, 32197806,
  "rs12903325",   15, 50061080,
  "rs1317082",     3, 169779797,
  "rs112540634",   6, 34656128,
  "rs7800079",     7, 139041601,
  "rs11596235",   10, 102631277,
  "rs10506957",   12, 88573230,
  "rs9535403",    13, 49898596,
  "rs7170049",    15, 41192089,
  "rs116403919",  17, 7167543,
  "rs11668706",   19, 10757583,
  "rs4918404",    10, 94567810,
  "rs4709820",     6, 164042540,
  "rs62550973",    9, 88788250,
  "rs62227732",   22, 30357033)

### African ancestry:
### Load data:
gwas_filtered <- vroom("GCST90165263_buildGRCh38.tsv.gz",
                       col_select = c(chromosome, base_pair_location, effect_allele, other_allele,
                                      odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency),
                       progress = TRUE) %>%
  mutate(chromosome = as.integer(chromosome),
         base_pair_location = as.integer(base_pair_location)) %>%
  semi_join(target_snps, by = c("chromosome" = "chr", "base_pair_location" = "pos"))

### Add rsIDs:
gwas_annotated <- gwas_filtered %>%
  left_join(target_snps, by = c("chromosome" = "chr", "base_pair_location" = "pos")) %>%
  relocate(rsid)

### Create columns:
gwas_annotated$beta <- log(gwas_annotated$odds_ratio)
gwas_annotated$se <- (((log(gwas_annotated$ci_upper))-(log(gwas_annotated$odds_ratio)))/1.96)
gwas_annotated <- gwas_annotated %>% select(rsid, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, se, p_value)
gwas_annotated <- gwas_annotated %>%
  rename(chr_afr = chromosome,
         bp_afr = base_pair_location,
         ea_afr = effect_allele,
         oa_afr = other_allele,
         eaf_afr = effect_allele_frequency,
         beta_afr = beta,
         se_afr = se,
         p_afr = p_value)

### Save:
write_tsv(gwas_annotated, "overall_inc_kessler_afr_GCST90165263_novel_snps.csv")

### East Asian ancestry:
### Load data:
gwas_filtered <- vroom("GCST90165265_buildGRCh38.tsv.gz",
                       col_select = c(chromosome, base_pair_location, effect_allele, other_allele,
                                      odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency),
                       progress = TRUE) %>%
  mutate(chromosome = as.integer(chromosome),
         base_pair_location = as.integer(base_pair_location)) %>%
  semi_join(target_snps, by = c("chromosome" = "chr", "base_pair_location" = "pos"))

### Add rsIDs:
gwas_annotated <- gwas_filtered %>%
  left_join(target_snps, by = c("chromosome" = "chr", "base_pair_location" = "pos")) %>%
  relocate(rsid)

### Create columns:
gwas_annotated$beta <- log(gwas_annotated$odds_ratio)
gwas_annotated$se <- (((log(gwas_annotated$ci_upper))-(log(gwas_annotated$odds_ratio)))/1.96)
gwas_annotated <- gwas_annotated %>% select(rsid, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, se, p_value)
gwas_annotated <- gwas_annotated %>%
  rename(chr_eas = chromosome,
         bp_eas = base_pair_location,
         ea_eas = effect_allele,
         oa_eas = other_allele,
         eaf_eas = effect_allele_frequency,
         beta_eas = beta,
         se_eas = se,
         p_eas = p_value)

### Save:
write_tsv(gwas_annotated, "overall_inc_kessler_eas_GCST90165265_novel_snps.csv")

### South Asian ancestry:
### Load data:
gwas_filtered <- vroom("GCST90165269_buildGRCh38.tsv.gz",
                       col_select = c(chromosome, base_pair_location, effect_allele, other_allele,
                                      odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency),
                       progress = TRUE) %>%
  mutate(chromosome = as.integer(chromosome),
         base_pair_location = as.integer(base_pair_location)) %>%
  semi_join(target_snps, by = c("chromosome" = "chr", "base_pair_location" = "pos"))

### Add rsIDs:
gwas_annotated <- gwas_filtered %>%
  left_join(target_snps, by = c("chromosome" = "chr", "base_pair_location" = "pos")) %>%
  relocate(rsid)

### Create columns:
gwas_annotated$beta <- log(gwas_annotated$odds_ratio)
gwas_annotated$se <- (((log(gwas_annotated$ci_upper))-(log(gwas_annotated$odds_ratio)))/1.96)
gwas_annotated <- gwas_annotated %>% select(rsid, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, se, p_value)
gwas_annotated <- gwas_annotated %>%
  rename(chr_eas = chromosome,
         bp_eas = base_pair_location,
         ea_eas = effect_allele,
         oa_eas = other_allele,
         eaf_eas = effect_allele_frequency,
         beta_eas = beta,
         se_eas = se,
         p_eas = p_value)

### Save:
write_tsv(gwas_annotated, "overall_inc_kessler_sas_GCST90165269_novel_snps.csv")