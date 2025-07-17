### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
ch_kessler_fp <- vroom("ch_kessler_sffdr_fuma_leadsnps_all_loci.csv")

### Calculate distances between sffdr loci and nearest kessler loci:
kessler_loci <- ch_kessler_fp %>%
  filter(grepl("Kessler", CH_GWAS))
sffdr_loci <- ch_kessler_fp %>%
  rowwise() %>%
  mutate(min_distance = if_else(any(kessler_loci$chr == chr),
                                min(abs(kessler_loci$pos[kessler_loci$chr == chr] - pos), na.rm = TRUE),
                                Inf),
         overlap = if_else(any(kessler_loci$chr == chr & kessler_loci$pos == pos), "yes", "no"),
         potentially_novel = if_else(!is.na(min_distance) & min_distance > 500000 & overlap == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, potentially_novel, everything())
sum(sffdr_loci$potentially_novel == "TRUE")
write_csv(sffdr_loci, file = "ch_kessler_fuma_loci_vs_sffdr_fuma_leadsnps.csv")

### Identify potentially novel loci:
sffdr_loci2 <- sffdr_loci %>%
  filter(potentially_novel == TRUE)
write_csv(sffdr_loci2, file = "ch_sffdr_fuma_potentially_novel_leadsnps_vs_kessler.csv")

### Identify potentially novel loci by CH type:
overall_inc <- sffdr_loci2 %>%
  filter(CH_GWAS == "sfFDR [Overall (inc)]")
dim(overall_inc)
overall_exc <- sffdr_loci2 %>%
  filter(CH_GWAS == "sfFDR [Overall (exc)]")
dim(overall_exc)
dnmt3a <- sffdr_loci2 %>%
  filter(CH_GWAS == "sfFDR [DNMT3A]")
dim(dnmt3a)
tet2 <- sffdr_loci2 %>%
  filter(CH_GWAS == "sfFDR [TET2]")
dim(tet2)
asxl1 <- sffdr_loci2 %>%
  filter(CH_GWAS == "sfFDR [ASXL1]")
dim(asxl1)

### Identify gene-specific novel loci not found in overall_inc:
overall_exc_novel <- overall_exc %>%
  rowwise() %>%
  mutate(min_distance_overall_inc = if_else(any(overall_inc$chr == chr),
                                            min(abs(overall_inc$pos[overall_inc$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_inc = if_else(any(overall_inc$chr == chr & overall_inc$pos == pos), "yes", "no"),
         potentially_novel_overall_inc = if_else(!is.na(min_distance_overall_inc) & min_distance_overall_inc > 500000 & overlap_overall_inc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
overall_exc_novel2 <- overall_exc_novel %>%
  filter(potentially_novel_overall_inc == "TRUE")
overall_exc_novel2

dnmt3a_novel <- dnmt3a %>%
  rowwise() %>%
  mutate(min_distance_overall_inc = if_else(any(overall_inc$chr == chr),
                                            min(abs(overall_inc$pos[overall_inc$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_inc = if_else(any(overall_inc$chr == chr & overall_inc$pos == pos), "yes", "no"),
         potentially_novel_overall_inc = if_else(!is.na(min_distance_overall_inc) & min_distance_overall_inc > 500000 & overlap_overall_inc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
dnmt3a_novel2 <- dnmt3a_novel %>%
  filter(potentially_novel_overall_inc == "TRUE")
dnmt3a_novel2

tet2_novel <- tet2 %>%
  rowwise() %>%
  mutate(min_distance_overall_inc = if_else(any(overall_inc$chr == chr),
                                            min(abs(overall_inc$pos[overall_inc$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_inc = if_else(any(overall_inc$chr == chr & overall_inc$pos == pos), "yes", "no"),
         potentially_novel_overall_inc = if_else(!is.na(min_distance_overall_inc) & min_distance_overall_inc > 500000 & overlap_overall_inc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
tet2_novel2 <- tet2_novel %>%
  filter(potentially_novel_overall_inc == "TRUE")
tet2_novel2

asxl1_novel <- asxl1 %>%
  rowwise() %>%
  mutate(min_distance_overall_inc = if_else(any(overall_inc$chr == chr),
                                            min(abs(overall_inc$pos[overall_inc$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_inc = if_else(any(overall_inc$chr == chr & overall_inc$pos == pos), "yes", "no"),
         potentially_novel_overall_inc = if_else(!is.na(min_distance_overall_inc) & min_distance_overall_inc > 500000 & overlap_overall_inc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
asxl1_novel2 <- asxl1_novel %>%
  filter(potentially_novel_overall_inc == "TRUE")
asxl1_novel2

### Identify gene-specific novel loci not found in overall_exc:
dnmt3a_novel3 <- dnmt3a_novel2 %>%
  rowwise() %>%
  mutate(min_distance_overall_exc = if_else(any(overall_exc_novel2$chr == chr),
                                            min(abs(overall_exc_novel2$pos[overall_exc_novel2$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_exc  = if_else(any(overall_exc_novel2$chr == chr & overall_exc_novel2$pos == pos), "yes", "no"),
         potentially_novel_overall_exc = if_else(!is.na(min_distance_overall_exc) & min_distance_overall_exc  > 500000 & overlap_overall_exc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
dnmt3a_novel4 <- dnmt3a_novel3 %>%
  filter(potentially_novel_overall_exc == "TRUE")
dnmt3a_novel4

tet2_novel3 <- tet2_novel2 %>%
  rowwise() %>%
  mutate(min_distance_overall_exc = if_else(any(overall_exc_novel2$chr == chr),
                                            min(abs(overall_exc_novel2$pos[overall_exc_novel2$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_exc  = if_else(any(overall_exc_novel2$chr == chr & overall_exc_novel2$pos == pos), "yes", "no"),
         potentially_novel_overall_exc = if_else(!is.na(min_distance_overall_exc) & min_distance_overall_exc  > 500000 & overlap_overall_exc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
tet2_novel4 <- tet2_novel3 %>%
  filter(potentially_novel_overall_exc == "TRUE")
tet2_novel4

asxl1_novel3 <- asxl1_novel2 %>%
  rowwise() %>%
  mutate(min_distance_overall_exc = if_else(any(overall_exc_novel2$chr == chr),
                                            min(abs(overall_exc_novel2$pos[overall_exc_novel2$chr == chr] - pos), na.rm = TRUE),
                                            Inf),
         overlap_overall_exc  = if_else(any(overall_exc_novel2$chr == chr & overall_exc_novel2$pos == pos), "yes", "no"),
         potentially_novel_overall_exc = if_else(!is.na(min_distance_overall_exc) & min_distance_overall_exc  > 500000 & overlap_overall_exc == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, min_distance, min_distance_overall_inc, overlap_overall_inc, potentially_novel_overall_inc, everything())
asxl1_novel4 <- asxl1_novel3 %>%
  filter(potentially_novel_overall_exc == "TRUE")
asxl1_novel4

### Select columns:
overall_inc <- overall_inc %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene)
overall_exc_novel2 <- overall_exc_novel2 %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene)
dnmt3a_novel4 <- dnmt3a_novel4 %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene)
tet2_novel4 <- tet2_novel4 %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene)
asxl1_novel4 <- asxl1_novel4 %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene)

### Save:
write_csv(overall_inc, "overall_inc_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
write_csv(overall_exc_novel2, file = "overall_exc_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
write_csv(dnmt3a_novel4, "dnmt3a_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
write_csv(tet2_novel4, "tet2_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
write_csv(asxl1_novel4, "asxl1_sffdr_fuma_novel_leadsnps_vs_kessler.csv")