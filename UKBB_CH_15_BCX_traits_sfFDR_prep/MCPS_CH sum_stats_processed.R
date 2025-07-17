### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)
library(MungeSumstats)

### Load and format g1000:
g1000 <- vroom(file = "g1000_eur.bim", delim = "\t", col_names = F)
g1000 <- dplyr::select(g1000, -c(X3, X5, X6))
g1000 <- rename(g1000, rsid = X2)
g1000 <- unite(g1000, "SNP", c(X1, X4), remove = TRUE)

### Load mcps_overall ch data:
mcps_overall <- vroom("GCST90435341.tsv")

### Format mcps_overall:
mcps_overall <- filter(mcps_overall, effect_allele_frequency >= 0.005)
mcps_overall <- filter(mcps_overall, effect_allele_frequency <= 0.995)
mcps_overall <- rename(mcps_overall, c(CHR = chromosome, BP = base_pair_location))
mcps_overall <- unite(mcps_overall, "SNP", CHR:BP, remove = FALSE)
mcps_overall <- data.table(mcps_overall)
mcps_overall_lifted <- liftover(mcps_overall, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
mcps_overall_lifted <- as_tibble(mcps_overall_lifted)
mcps_overall_lifted <- select(mcps_overall_lifted, -SNP)
mcps_overall_lifted <- unite(mcps_overall_lifted, "SNP", CHR:BP, remove = FALSE)
mcps_overall_rsid <- inner_join(mcps_overall_lifted, g1000, by = "SNP")
mcps_overall_rsid <- mcps_overall_rsid %>% mutate(p_value = 10^(-neg_log_10_p_value))
mcps_overall_rsid <- mcps_overall_rsid %>% distinct()
mcps_overall_rsid <- mcps_overall_rsid %>% select(SNP, rsid, CHR, BP, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)
mcps_overall_rsid <- mcps_overall_rsid %>%
  rename(ea_mcps_overall = effect_allele,
         oa_mcps_overall = other_allele,
         eaf_mcps_overall = effect_allele_frequency,
         beta_mcps_overall = beta,
         se_mcps_overall = standard_error,
         pval_mcps_overall = p_value)
vroom_write(mcps_overall_rsid, file = "mcps_overall_GCST90435341_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load mcps_dnmt3a:
mcps_dnmt3a <- vroom("GCST90435342.tsv")

### Format mcps_dnmt3a:
mcps_dnmt3a <- filter(mcps_dnmt3a, effect_allele_frequency >= 0.005)
mcps_dnmt3a <- filter(mcps_dnmt3a, effect_allele_frequency <= 0.995)
mcps_dnmt3a <- rename(mcps_dnmt3a, c(CHR = chromosome, BP = base_pair_location))
mcps_dnmt3a <- unite(mcps_dnmt3a, "SNP", CHR:BP, remove = FALSE)
mcps_dnmt3a <- data.table(mcps_dnmt3a)
mcps_dnmt3a_lifted <- liftover(mcps_dnmt3a, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
mcps_dnmt3a_lifted <- as_tibble(mcps_dnmt3a_lifted)
mcps_dnmt3a_lifted <- select(mcps_dnmt3a_lifted, -SNP)
mcps_dnmt3a_lifted <- unite(mcps_dnmt3a_lifted, "SNP", CHR:BP, remove = FALSE)
mcps_dnmt3a_rsid <- inner_join(mcps_dnmt3a_lifted, g1000, by = "SNP")
mcps_dnmt3a_rsid <- mcps_dnmt3a_rsid %>% mutate(p_value = 10^(-neg_log_10_p_value))
mcps_dnmt3a_rsid <- mcps_dnmt3a_rsid %>% distinct()
mcps_dnmt3a_rsid <- mcps_dnmt3a_rsid %>% select(SNP, rsid, CHR, BP, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)
mcps_dnmt3a_rsid <- mcps_dnmt3a_rsid %>%
  rename(ea_mcps_dnmt3a = effect_allele,
         oa_mcps_dnmt3a = other_allele,
         eaf_mcps_dnmt3a = effect_allele_frequency,
         beta_mcps_dnmt3a = beta,
         se_mcps_dnmt3a = standard_error,
         pval_mcps_dnmt3a = p_value)
vroom_write(mcps_dnmt3a_rsid, file = "mcps_dnmt3a_GCST90435342_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load mcps_tet2:
mcps_tet2 <- vroom("GCST90435343.tsv")

### Format mcps_tet2:
mcps_tet2 <- filter(mcps_tet2, effect_allele_frequency >= 0.005)
mcps_tet2 <- filter(mcps_tet2, effect_allele_frequency <= 0.995)
mcps_tet2 <- rename(mcps_tet2, c(CHR = chromosome, BP = base_pair_location))
mcps_tet2 <- unite(mcps_tet2, "SNP", CHR:BP, remove = FALSE)
mcps_tet2 <- data.table(mcps_tet2)
mcps_tet2_lifted <- liftover(mcps_tet2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
mcps_tet2_lifted <- as_tibble(mcps_tet2_lifted)
mcps_tet2_lifted <- select(mcps_tet2_lifted, -SNP)
mcps_tet2_lifted <- unite(mcps_tet2_lifted, "SNP", CHR:BP, remove = FALSE)
mcps_tet2_rsid <- inner_join(mcps_tet2_lifted, g1000, by = "SNP")
mcps_tet2_rsid <- mcps_tet2_rsid %>% mutate(p_value = 10^(-neg_log_10_p_value))
mcps_tet2_rsid <- mcps_tet2_rsid %>% distinct()
mcps_tet2_rsid <- mcps_tet2_rsid %>% select(SNP, rsid, CHR, BP, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)
mcps_tet2_rsid <- mcps_tet2_rsid %>%
  rename(ea_mcps_tet2 = effect_allele,
         oa_mcps_tet2 = other_allele,
         eaf_mcps_tet2 = effect_allele_frequency,
         beta_mcps_tet2 = beta,
         se_mcps_tet2 = standard_error,
         pval_mcps_tet2 = p_value)
vroom_write(mcps_tet2_rsid, file = "mcps_tet2_GCST90435343_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Load mcps_asxl1:
mcps_asxl1 <- vroom("GCST90435344.tsv")

### Format mcps_asxl1:
mcps_asxl1 <- filter(mcps_asxl1, effect_allele_frequency >= 0.005)
mcps_asxl1 <- filter(mcps_asxl1, effect_allele_frequency <= 0.995)
mcps_asxl1 <- rename(mcps_asxl1, c(CHR = chromosome, BP = base_pair_location))
mcps_asxl1 <- unite(mcps_asxl1, "SNP", CHR:BP, remove = FALSE)
mcps_asxl1 <- data.table(mcps_asxl1)
mcps_asxl1_lifted <- liftover(mcps_asxl1, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
mcps_asxl1_lifted <- as_tibble(mcps_asxl1_lifted)
mcps_asxl1_lifted <- select(mcps_asxl1_lifted, -SNP)
mcps_asxl1_lifted <- unite(mcps_asxl1_lifted, "SNP", CHR:BP, remove = FALSE)
mcps_asxl1_rsid <- inner_join(mcps_asxl1_lifted, g1000, by = "SNP")
mcps_asxl1_rsid <- mcps_asxl1_rsid %>% mutate(p_value = 10^(-neg_log_10_p_value))
mcps_asxl1_rsid <- mcps_asxl1_rsid %>% distinct()
mcps_asxl1_rsid <- mcps_asxl1_rsid %>% select(SNP, rsid, CHR, BP, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)
mcps_asxl1_rsid <- mcps_asxl1_rsid %>%
  rename(ea_mcps_asxl1 = effect_allele,
         oa_mcps_asxl1 = other_allele,
         eaf_mcps_asxl1 = effect_allele_frequency,
         beta_mcps_asxl1 = beta,
         se_mcps_asxl1 = standard_error,
         pval_mcps_asxl1 = p_value)
vroom_write(mcps_asxl1_rsid, file = "mcps_asxl1_GCST90435344_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")