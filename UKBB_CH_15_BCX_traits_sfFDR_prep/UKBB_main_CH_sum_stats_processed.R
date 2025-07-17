### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)
library(MungeSumstats)

### Load data:
overall <- vroom("GCST90165261_buildGRCh38.tsv.gz",  col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))
g1000 <- vroom(file = "g1000_eur.bim", delim = "\t", col_names = F)

### Format g1000:
g1000 <- dplyr::select(g1000, -c(X3, X5, X6))
g1000 <- rename(g1000, rsid = X2)
g1000 <- unite(g1000, "SNP", c(X1, X4), remove = TRUE)

### Format overall:
overall <- filter(overall, effect_allele_frequency >= 0.005)
overall <- filter(overall, effect_allele_frequency <= 0.995)
overall <- rename(overall, c(CHR = chromosome, BP = base_pair_location))
overall <- unite(overall, "SNP", CHR:BP, remove = FALSE)
overall <- data.table(overall)
overall_lifted <- liftover(overall, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
overall_lifted <- as_tibble(overall_lifted)
overall_lifted <- select(overall_lifted, -SNP)
overall_lifted <- unite(overall_lifted, "SNP", CHR:BP, remove = FALSE)
overall_rsid <- inner_join(overall_lifted, g1000, by = "SNP")
overall_rsid$beta <- log(overall_rsid$odds_ratio)
overall_rsid$se <- (((log(overall_rsid$ci_upper))-(log(overall_rsid$odds_ratio)))/1.96)
overall_rsid <- overall_rsid %>%
 distinct()
overall_rsid <- overall_rsid %>%
 select(SNP, CHR, BP, other_allele, effect_allele, p_value, effect_allele_frequency, rsid, beta, se)
overall_rsid <- overall_rsid %>%
  rename(ea_overall = effect_allele,
         oa_overall = other_allele,
         eaf_overall = effect_allele_frequency,
         beta_overall = beta,
         se_overall = se,
         pval_overall = p_value)
vroom_write(overall_rsid, file = "overall_GCST90165261_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load data:
overall <- vroom("GCST90165267_buildGRCh38.tsv.gz",  col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))

### Format overall:
overall <- filter(overall, effect_allele_frequency >= 0.005)
overall <- filter(overall, effect_allele_frequency <= 0.995)
overall <- rename(overall, c(CHR = chromosome, BP = base_pair_location))
overall <- unite(overall, "SNP", CHR:BP, remove = FALSE)
overall <- data.table(overall)
overall_lifted <- liftover(overall, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
overall_lifted <- as_tibble(overall_lifted)
overall_lifted <- select(overall_lifted, -SNP)
overall_lifted <- unite(overall_lifted, "SNP", CHR:BP, remove = FALSE)
overall_rsid <- inner_join(overall_lifted, g1000, by = "SNP")
overall_rsid$beta <- log(overall_rsid$odds_ratio)
overall_rsid$se <- (((log(overall_rsid$ci_upper))-(log(overall_rsid$odds_ratio)))/1.96)
overall_rsid <- overall_rsid %>% distinct()
overall_rsid <- overall_rsid %>% select(SNP, CHR, BP, other_allele, effect_allele, p_value, effect_allele_frequency, rsid, beta, se)
overall_rsid <- overall_rsid %>%
  rename(ea_overall = effect_allele,
         oa_overall = other_allele,
         eaf_overall = effect_allele_frequency,
         beta_overall = beta,
         se_overall = se,
         pval_overall = p_value)
vroom_write(overall_rsid, file = "overall_GCST90165267_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load dnmt3a data:
dnmt3a <- vroom("GCST90165271_buildGRCh38.tsv.gz",  col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))

### Format dnmt3a:
dnmt3a <- filter(dnmt3a, effect_allele_frequency >= 0.005)
dnmt3a <- filter(dnmt3a, effect_allele_frequency <= 0.995)
dnmt3a <- rename(dnmt3a, c(CHR = chromosome, BP = base_pair_location))
dnmt3a <- unite(dnmt3a, "SNP", CHR:BP, remove = FALSE)
dnmt3a <- data.table(dnmt3a)
dnmt3a_lifted <- liftover(dnmt3a, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
dnmt3a_lifted <- as_tibble(dnmt3a_lifted)
dnmt3a_lifted <- select(dnmt3a_lifted, -SNP)
dnmt3a_lifted <- unite(dnmt3a_lifted, "SNP", CHR:BP, remove = FALSE)
dnmt3a_rsid <- inner_join(dnmt3a_lifted, g1000, by = "SNP")
dnmt3a_rsid$beta <- log(dnmt3a_rsid$odds_ratio)
dnmt3a_rsid$se <- (((log(dnmt3a_rsid$ci_upper))-(log(dnmt3a_rsid$odds_ratio)))/1.96)
dnmt3a_rsid <- dnmt3a_rsid %>% distinct()
dnmt3a_rsid <- dnmt3a_rsid %>% select(SNP, CHR, BP, other_allele, effect_allele, p_value, effect_allele_frequency, rsid, beta, se)
dnmt3a_rsid <- dnmt3a_rsid %>%
  rename(ea_dnmt3a = effect_allele,
         oa_dnmt3a = other_allele,
         eaf_dnmt3a = effect_allele_frequency,
         beta_dnmt3a = beta,
         se_dnmt3a = se,
         pval_dnmt3a = p_value)
vroom_write(dnmt3a_rsid, file = "dnmt3a_GCST90165271_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load tet2 data:
tet2 <- vroom("GCST90165281_buildGRCh38.tsv.gz", col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))

### Format tet2:
tet2 <- filter(tet2, effect_allele_frequency >= 0.005)
tet2 <- filter(tet2, effect_allele_frequency <= 0.995)
tet2 <- rename(tet2, c(CHR = chromosome, BP = base_pair_location))
tet2 <- unite(tet2, "SNP", CHR:BP, remove = FALSE)
tet2 <- data.table(tet2)
tet2_lifted <- liftover(tet2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
tet2_lifted <- as_tibble(tet2_lifted)
tet2_lifted <- select(tet2_lifted, -SNP)
tet2_lifted <- unite(tet2_lifted, "SNP", CHR:BP, remove = FALSE)
tet2_rsid <- inner_join(tet2_lifted, g1000, by = "SNP")
tet2_rsid$beta <- log(tet2_rsid$odds_ratio)
tet2_rsid$se <- (((log(tet2_rsid$ci_upper))-(log(tet2_rsid$odds_ratio)))/1.96)
tet2_rsid <- tet2_rsid %>% distinct()
tet2_rsid <- tet2_rsid %>% select(SNP, CHR, BP, other_allele, effect_allele, p_value, effect_allele_frequency, rsid, beta, se)
tet2_rsid <- tet2_rsid %>%
  rename(ea_tet2 = effect_allele,
         oa_tet2 = other_allele,
         eaf_tet2 = effect_allele_frequency,
         beta_tet2 = beta,
         se_tet2 = se,
         pval_tet2 = p_value)
vroom_write(tet2_rsid, file = "tet2_GCST90165281_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load asxl1 data:
asxl1 <- vroom("GCST90165259_buildGRCh38.tsv.gz" , col_select = c(chromosome, base_pair_location, other_allele, effect_allele, odds_ratio, ci_lower, ci_upper, p_value, effect_allele_frequency))

### Format asxl1:
asxl1 <- filter(asxl1, effect_allele_frequency >= 0.005)
asxl1 <- filter(asxl1, effect_allele_frequency <= 0.995)
asxl1 <- rename(asxl1, c(CHR = chromosome, BP = base_pair_location))
asxl1 <- unite(asxl1, "SNP", CHR:BP, remove = FALSE)
asxl1 <- data.table(asxl1)
asxl1_lifted <- liftover(asxl1, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
asxl1_lifted <- as_tibble(asxl1_lifted)
asxl1_lifted <- select(asxl1_lifted, -SNP)
asxl1_lifted <- unite(asxl1_lifted, "SNP", CHR:BP, remove = FALSE)
asxl1_rsid <- inner_join(asxl1_lifted, g1000, by = "SNP")
asxl1_rsid$beta <- log(asxl1_rsid$odds_ratio)
asxl1_rsid$se <- (((log(asxl1_rsid$ci_upper))-(log(asxl1_rsid$odds_ratio)))/1.96)
asxl1_rsid <- asxl1_rsid %>% distinct()
asxl1_rsid <- asxl1_rsid %>% select(SNP, CHR, BP, other_allele, effect_allele, p_value, effect_allele_frequency, rsid, beta, se)
asxl1_rsid <- asxl1_rsid %>%
  rename(ea_asxl1 = effect_allele,
         oa_asxl1 = other_allele,
         eaf_asxl1 = effect_allele_frequency,
         beta_asxl1 = beta,
         se_asxl1 = se,
         pval_asxl1 = p_value)
vroom_write(asxl1_rsid, file = "asxl1_GCST90165259_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")