### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load data:
overall_inc_kessler <- vroom("overall_inc_kessler_fuma_leadsnps_nearest_gene.csv")
overall_exc_kessler <- vroom("overall_exc_kessler_fuma_leadsnps_nearest_gene.csv")
dnmt3a_kessler <- vroom("dnmt3a_kessler_fuma_leadsnps_nearest_gene.csv")
tet2_kessler <- vroom("tet2_kessler_fuma_leadsnps_nearest_gene.csv")
asxl1_kessler <- vroom("asxl1_kessler_fuma_leadsnps_nearest_gene.csv")

overall_inc_fp <- vroom("overall_inc_fp_fuma_leadsnps_nearest_gene.csv")
overall_exc_fp <- vroom("overall_exc_fp_fuma_leadsnps_nearest_gene.csv")
dnmt3a_fp <- vroom("dnmt3a_fp_fuma_leadsnps_nearest_gene.csv")
tet2_fp <- vroom("tet2_fp_fuma_leadsnps_nearest_gene.csv")
asxl1_fp <- vroom("asxl1_fp_fuma_leadsnps_nearest_gene.csv")

### Format kessler data:
overall_inc_kessler2 <- overall_inc_kessler %>%
  mutate(CH_GWAS = "Kessler [Overall (inc)]")
overall_exc_kessler2 <- overall_exc_kessler %>%
  mutate(CH_GWAS = "Kessler [Overall (exc)]")
dnmt3a_kessler2 <- dnmt3a_kessler %>%
  mutate(CH_GWAS = "Kessler [DNMT3A]")
tet2_kessler2 <- tet2_kessler %>%
  mutate(CH_GWAS = "Kessler [TET2]]")
asxl1_kessler2 <- asxl1_kessler %>%
  mutate(CH_GWAS = "Kessler [ASXL1]")

### Format sffdr data:
overall_inc_fp2 <- overall_inc_fp %>%
  mutate(CH_GWAS = "sfFDR [Overall (inc)]")
overall_exc_fp2 <- overall_exc_fp %>%
  mutate(CH_GWAS = "sfFDR [Overall (exc)]")
dnmt3a_fp2 <- dnmt3a_fp %>%
  mutate(CH_GWAS = "sfFDR [DNMT3A]")
tet2_fp2 <- tet2_fp %>%
  mutate(CH_GWAS = "sfFDR [TET2]")
asxl1_fp2 <- asxl1_fp %>%
  mutate(CH_GWAS = "sfFDR [ASXL1]")

### Combine all data:
ch_kessler_fp <- rbind(overall_inc_kessler2,
                       overall_exc_kessler2, 
                       dnmt3a_kessler2, 
                       tet2_kessler2, 
                       asxl1_kessler2, 
                       overall_inc_fp2, 
                       overall_exc_fp2, 
                       dnmt3a_fp2, 
                       tet2_fp2, 
                       asxl1_fp2)
dim(ch_kessler_fp)

### Format ch_kessler_fp:
ch_kessler_fp2 <- ch_kessler_fp %>%
  arrange(chr, pos) %>%
  select(CH_GWAS, GenomicLocus, rsID, chr, pos, uniqID, nIndSigSNPs, IndSigSNPs, p, nearestGene)

### Save output:
write_csv(ch_kessler_fp2, "ch_kessler_sffdr_fuma_leadsnps_all_loci.csv")

### Filter all sffdr loci:
ch_kessler_fp3 <-  ch_kessler_fp2 %>%
 filter(str_detect(CH_GWAS, "sfFDR"))

### Save output:
write_csv(ch_kessler_fp3, "sffdr_all_ch_gws_snps.csv")
