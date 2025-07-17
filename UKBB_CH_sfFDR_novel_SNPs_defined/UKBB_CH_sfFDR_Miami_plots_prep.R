### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_inc_kessler <- vroom("overall_inc_kessler_orig_fuma.txt.gz")
overall_exc_kessler <- vroom("overall_exc_kessler_orig_fuma.txt.gz")
dnmt3a_kessler <- vroom("dnmt3a_kessler_orig_fuma.txt.gz")
tet2_kessler <- vroom("tet2_kessler_orig_fuma.txt.gz")
asxl1_kessler <- vroom("asxl1_kessler_orig_fuma.txt.gz")

overall_inc_kessler_genes <- vroom("overall_inc_kessler_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
overall_exc_kessler_genes <- vroom("overall_exc_kessler_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
dnmt3a_kessler_genes <- vroom("dnmt3a_kessler_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
tet2_kessler_genes <- vroom("tet2_kessler_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
asxl1_kessler_genes <- vroom("asxl1_kessler_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))

overall_inc_sffdr <- vroom("overall_inc_kessler_sffdr_fp_fuma.txt.gz")
overall_exc_sffdr <- vroom("overall_exc_kessler_sffdr_fp_fuma.txt.gz")
dnmt3a_sffdr <- vroom("dnmt3a_kessler_sffdr_fp_fuma.txt.gz")
tet2_sffdr <- vroom("tet2_kessler_sffdr_fp_fuma.txt.gz")
asxl1_sffdr <- vroom("asxl1_kessler_sffdr_fp_fuma.txt.gz")

overall_inc_sffdr_genes <- vroom("overall_inc_fp_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
overall_exc_sffdr_genes <- vroom("overall_exc_fp_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
dnmt3a_sffdr_genes <- vroom("dnmt3a_fp_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
tet2_sffdr_genes <- vroom("tet2_fp_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))
asxl1_sffdr_genes <- vroom("asxl1_fp_fuma_leadsnps_nearest_gene.csv", col_select = c(rsID, nearestGene))

### Rename SNP:
overall_inc_kessler_genes <- overall_inc_kessler_genes %>%
  rename(SNP = rsID)
overall_exc_kessler_genes <- overall_exc_kessler_genes %>%
  rename(SNP = rsID)
dnmt3a_kessler_genes <- dnmt3a_kessler_genes %>%
  rename(SNP = rsID)
tet2_kessler_genes <- tet2_kessler_genes %>%
  rename(SNP = rsID)
asxl1_kessler_genes <- asxl1_kessler_genes %>%
  rename(SNP = rsID)
overall_inc_sffdr_genes <- overall_inc_sffdr_genes %>%
  rename(SNP = rsID)
overall_exc_sffdr_genes <- overall_exc_sffdr_genes %>%
  rename(SNP = rsID)
dnmt3a_sffdr_genes <- dnmt3a_sffdr_genes %>%
  rename(SNP = rsID)
tet2_sffdr_genes <- tet2_sffdr_genes %>%
  rename(SNP = rsID)
asxl1_sffdr_genes <- asxl1_sffdr_genes %>%
  rename(SNP = rsID)

### Format original data for plotting:
overall_inc_kessler2 <- overall_inc_kessler %>%
  left_join(overall_inc_kessler_genes, by = "SNP") %>%
  mutate(study = "overall_inc_kessler") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

overall_exc_kessler2 <- overall_exc_kessler %>%
  left_join(overall_exc_kessler_genes, by = "SNP") %>%
  mutate(study = "overall_exc_kessler") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

dnmt3a_kessler2 <- dnmt3a_kessler %>%
  left_join(dnmt3a_kessler_genes, by = "SNP") %>%
  mutate(study = "dnmt3a_kessler") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

tet2_kessler2 <- tet2_kessler %>%
  left_join(tet2_kessler_genes, by = "SNP") %>%
  mutate(study = "tet2_kessler") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

asxl1_kessler2 <- asxl1_kessler %>%
  left_join(overall_inc_kessler_genes, by = "SNP") %>%
  mutate(study = "asxl1_kessler") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

### Format sffdr data for plotting:
overall_inc_sffdr2 <- overall_inc_sffdr %>%
  left_join(overall_inc_sffdr_genes, by = "SNP") %>%
  mutate(study = "overall_inc_sffdr") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

overall_exc_sffdr2 <- overall_exc_sffdr %>%
  left_join(overall_exc_sffdr_genes, by = "SNP") %>%
  mutate(study = "overall_exc_sffdr") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

dnmt3a_sffdr2 <- dnmt3a_sffdr %>%
  left_join(dnmt3a_sffdr_genes, by = "SNP") %>%
  mutate(study = "dnmt3a_sffdr") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

tet2_sffdr2 <- tet2_sffdr %>%
  left_join(tet2_sffdr_genes, by = "SNP") %>%
  mutate(study = "tet2_sffdr") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

asxl1_sffdr2 <- asxl1_sffdr %>%
  left_join(asxl1_sffdr_genes, by = "SNP") %>%
  mutate(study = "asxl1_sffdr") %>%
  mutate(tstat = Beta / SE) %>%
  rename(rsid = SNP,
         chr = CHR,
         pos = BP,
         beta = Beta,
         se = SE,
         pval = P) %>%
  select(rsid, chr, pos, beta, se, tstat, pval, study, nearestGene) %>%
  filter(pval<0.05)

### Combine by ch subtype:
overall_inc_sffdr_kessler <- as.data.frame(bind_rows(overall_inc_sffdr2, overall_inc_kessler2))
overall_exc_sffdr_kessler <- as.data.frame(bind_rows(overall_exc_sffdr2, overall_exc_kessler2))
dnmt3a_sffdr_kessler <- as.data.frame(bind_rows(dnmt3a_sffdr2, dnmt3a_kessler2))
tet2_sffdr_kessler <- as.data.frame(bind_rows(tet2_sffdr2, tet2_kessler2))
asxl1_sffdr_kessler <- as.data.frame(bind_rows(asxl1_sffdr2, asxl1_kessler2))

### Save:
write_csv(overall_inc_sffdr_kessler, file = "overall_inc_sffdr_kessler_for_miami.csv" )
write_csv(overall_exc_sffdr_kessler, file = "overall_exc_sffdr_kessler_for_miami.csv")
write_csv(dnmt3a_sffdr_kessler, file = "dnmt3a_sffdr_kessler_for_miami.csv")
write_csv(tet2_sffdr_kessler, file = "tet2_sffdr_kessler_for_miami.csv")
write_csv(asxl1_sffdr_kessler, file = "asxl1_sffdr_kessler_for_miami.csv")  