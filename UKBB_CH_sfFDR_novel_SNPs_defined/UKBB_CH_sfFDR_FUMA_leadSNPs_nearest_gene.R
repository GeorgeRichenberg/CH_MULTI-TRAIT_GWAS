### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
overall_inc_kessler <- vroom("overall_inc_kessler_leadSNPs.txt")
overall_exc_kessler <- vroom("overall_exc_kessler_leadSNPs.txt")
dnmt3a_kessler <- vroom("dnmt3a_kessler_leadSNPs.txt")
tet2_kessler <- vroom("tet2_kessler_leadSNPs.txt")
asxl1_kessler <- vroom("asxl1_kessler_leadSNPs.txt")

overall_inc_fp <- vroom("overall_inc_fp_leadSNPs.txt")
overall_exc_fp <- vroom("overall_exc_fp_leadSNPs.txt")
dnmt3a_fp <- vroom("dnmt3a_fp_leadSNPs.txt")
tet2_fp <- vroom("tet2_fp_leadSNPs.txt")
asxl1_fp <- vroom("asxl1_fp_leadSNPs.txt")

overall_inc_kessler_genes <- vroom("overall_inc_kessler_snps.txt")
overall_exc_kessler_genes <- vroom("overall_exc_kessler_snps.txt")
dnmt3a_kessler_genes <- vroom("dnmt3a_kessler_snps.txt")
tet2_kessler_genes <- vroom("tet2_kessler_snps.txt")
asxl1_kessler_genes <- vroom("asxl1_kessler_snps.txt")

overall_inc_fp_genes <- vroom("overall_inc_fp_snps.txt")
overall_exc_fp_genes <- vroom("overall_exc_fp_snps.txt")
dnmt3a_fp_genes <- vroom("dnmt3a_fp_snps.txt")
tet2_fp_genes <- vroom("tet2_fp_snps.txt")
asxl1_fp_genes <- vroom("asxl1_fp_snps.txt")

### Format genes:
overall_inc_kessler_genes2 <- overall_inc_kessler_genes %>%
  select(uniqID, nearestGene)
overall_exc_kessler_genes2 <- overall_exc_kessler_genes %>%
  select(uniqID, nearestGene)
dnmt3a_kessler_genes2 <- dnmt3a_kessler_genes %>%
  select(uniqID, nearestGene)
tet2_kessler_genes2 <- tet2_kessler_genes %>%
  select(uniqID, nearestGene)
asxl1_kessler_genes2 <- asxl1_kessler_genes %>%
  select(uniqID, nearestGene)

overall_inc_fp_genes2 <- overall_inc_fp_genes %>%
  select(uniqID, nearestGene)
overall_exc_fp_genes2 <- overall_exc_fp_genes %>%
  select(uniqID, nearestGene)
dnmt3a_fp_genes2 <- dnmt3a_fp_genes %>%
  select(uniqID, nearestGene)
tet2_fp_genes2 <- tet2_fp_genes %>%
  select(uniqID, nearestGene)
asxl1_fp_genes2 <- asxl1_fp_genes %>%
  select(uniqID, nearestGene)

### Identify nearest gene for lead SNP:
overall_inc_kessler2 <- overall_inc_kessler %>%
  left_join(overall_inc_kessler_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(overall_inc_kessler)
dim(overall_inc_kessler2)

overall_exc_kessler2 <- overall_exc_kessler %>%
  left_join(overall_exc_kessler_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(overall_exc_kessler)
dim(overall_exc_kessler2)

dnmt3a_kessler2 <- dnmt3a_kessler %>%
  left_join(dnmt3a_kessler_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(dnmt3a_kessler)
dim(dnmt3a_kessler2)

tet2_kessler2 <- tet2_kessler %>%
  left_join(tet2_kessler_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(tet2_kessler)
dim(tet2_kessler2)

asxl1_kessler2 <- asxl1_kessler %>%
  left_join(asxl1_kessler_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(asxl1_kessler)
dim(asxl1_kessler2)

overall_inc_fp2 <- overall_inc_fp %>%
  left_join(overall_inc_fp_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(overall_inc_fp)
dim(overall_inc_fp2)

overall_exc_fp2 <- overall_exc_fp %>%
  left_join(overall_exc_fp_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(overall_exc_fp)
dim(overall_exc_fp2)

dnmt3a_fp2 <- dnmt3a_fp %>%
  left_join(dnmt3a_fp_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(dnmt3a_fp)
dim(dnmt3a_fp2)

tet2_fp2 <- tet2_fp %>%
  left_join(tet2_fp_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(tet2_fp)
dim(tet2_fp2)

asxl1_fp2 <- asxl1_fp %>%
  left_join(asxl1_fp_genes2, by = "uniqID") %>%
  mutate(nearestGene = if_else(is.na(nearestGene), "NA", nearestGene))
dim(asxl1_fp)
dim(asxl1_fp2)

### Save outputs:
write_csv(overall_inc_kessler2, file = "overall_inc_kessler_fuma_leadsnps_nearest_gene.csv")
write_csv(overall_exc_kessler2, file = "overall_exc_kessler_fuma_leadsnps_nearest_gene.csv")
write_csv(dnmt3a_kessler2, file = "dnmt3a_kessler_fuma_leadsnps_nearest_gene.csv")
write_csv(tet2_kessler2, file = "tet2_kessler_fuma_leadsnps_nearest_gene.csv")
write_csv(asxl1_kessler2, file = "asxl1_kessler_fuma_leadsnps_nearest_gene.csv")

write_csv(overall_inc_fp2, file = "overall_inc_fp_fuma_leadsnps_nearest_gene.csv")
write_csv(overall_exc_fp2, file = "overall_exc_fp_fuma_leadsnps_nearest_gene.csv")
write_csv(dnmt3a_fp2, file = "dnmt3a_fp_fuma_leadsnps_nearest_gene.csv")
write_csv(tet2_fp2, file = "tet2_fp_fuma_leadsnps_nearest_gene.csv")
write_csv(asxl1_fp2, file = "asxl1_fp_fuma_leadsnps_nearest_gene.csv")
