### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)
library(MungeSumstats)

### Load data:
wen_mcps <- vroom("wen_et_al_mcps_lead_snps.csv")
wen_meta <- vroom("wen_et_al_ch_gwas_lead_snps.csv")
overall_inc_novel <- vroom("overall_inc_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
overall_exc_novel <- vroom("overall_exc_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
dnmt3a_novel <- vroom("dnmt3a_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
tet2_novel <- vroom("tet2_sffdr_fuma_novel_leadsnps_vs_kessler.csv")
asxl1_novel <- vroom("asxl1_sffdr_fuma_novel_leadsnps_vs_kessler.csv")

### Separate wen_mcps by CH type:
wen_mcps_dnmt3a <- wen_mcps %>%
  filter(CH %in% c("All CH", "DNMT3A"))
wen_mcps_tet2 <- wen_mcps %>%
  filter(CH %in% c("All CH", "TET2"))
wen_mcps_asxl1 <- wen_mcps %>%
  filter(CH %in% c("All CH", "ASXL1"))

### Identify potentially novel loci overlapping with wen_mcps:
overall_inc_novel2 <- overall_inc_novel %>%
  rowwise() %>%
  mutate(min_distance_wen_mcps = if_else(any(wen_mcps$Chr == chr),
                                         min(abs(wen_mcps$Pos[wen_mcps$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_mcps = if_else(any(wen_mcps$Chr == chr & wen_mcps$Pos == pos), "yes", "no"),
         potentially_novel_wen_mcps = if_else(!is.na(min_distance_wen_mcps) & min_distance_wen_mcps > 500000 & overlap_wen_mcps == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, potentially_novel, min_distance_wen_mcps, overlap_wen_mcps, potentially_novel_wen_mcps, everything())
sum(overall_inc_novel2$potentially_novel_wen_mcps == "TRUE")

overall_exc_novel2 <- overall_exc_novel %>%
  rowwise() %>%
  mutate(min_distance_wen_mcps = if_else(any(wen_mcps$Chr == chr),
                                         min(abs(wen_mcps$Pos[wen_mcps$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_mcps = if_else(any(wen_mcps$Chr == chr & wen_mcps$Pos == pos), "yes", "no"),
         potentially_novel_wen_mcps = if_else(!is.na(min_distance_wen_mcps) & min_distance_wen_mcps > 500000 & overlap_wen_mcps == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, potentially_novel, min_distance_wen_mcps, overlap_wen_mcps, potentially_novel_wen_mcps, everything())
sum(overall_exc_novel2$potentially_novel_wen_mcps == "TRUE")

dnmt3a_novel2 <- dnmt3a_novel %>%
  rowwise() %>%
  mutate(min_distance_wen_mcps = if_else(any(wen_mcps_dnmt3a$Chr == chr),
                                         min(abs(wen_mcps_dnmt3a$Pos[wen_mcps_dnmt3a$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_mcps = if_else(any(wen_mcps_dnmt3a$Chr == chr & wen_mcps_dnmt3a$Pos == pos), "yes", "no"),
         potentially_novel_wen_mcps = if_else(!is.na(min_distance_wen_mcps) & min_distance_wen_mcps > 500000 & overlap_wen_mcps == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, potentially_novel, min_distance_wen_mcps, overlap_wen_mcps, potentially_novel_wen_mcps, everything())
sum(dnmt3a_novel2$potentially_novel_wen_mcps == "TRUE")

tet2_novel2 <- tet2_novel %>%
  rowwise() %>%
  mutate(min_distance_wen_mcps = if_else(any(wen_mcps_tet2$Chr == chr),
                                         min(abs(wen_mcps_tet2$Pos[wen_mcps_tet2$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_mcps = if_else(any(wen_mcps_tet2$Chr == chr & wen_mcps_tet2$Pos == pos), "yes", "no"),
         potentially_novel_wen_mcps = if_else(!is.na(min_distance_wen_mcps) & min_distance_wen_mcps > 500000 & overlap_wen_mcps == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, potentially_novel, min_distance_wen_mcps, overlap_wen_mcps, potentially_novel_wen_mcps, everything())
sum(tet2_novel2$potentially_novel_wen_mcps == "TRUE")

asxl1_novel2 <- asxl1_novel %>%
  rowwise() %>%
  mutate(min_distance_wen_mcps = if_else(any(wen_mcps_asxl1$Chr == chr),
                                         min(abs(wen_mcps_asxl1$Pos[wen_mcps_asxl1$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_mcps = if_else(any(wen_mcps_asxl1$Chr == chr & wen_mcps_asxl1$Pos == pos), "yes", "no"),
         potentially_novel_wen_mcps = if_else(!is.na(min_distance_wen_mcps) & min_distance_wen_mcps > 500000 & overlap_wen_mcps == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, potentially_novel, min_distance_wen_mcps, overlap_wen_mcps, potentially_novel_wen_mcps, everything())
sum(asxl1_novel2$potentially_novel_wen_mcps == "TRUE")

### Filter potentially_novel_wen_mcps = TRUE:
overall_inc_novel2 <- overall_inc_novel2 %>%
  filter(potentially_novel_wen_mcps == "TRUE")
dim(overall_inc_novel)
dim(overall_inc_novel2)
overall_exc_novel2 <- overall_exc_novel2 %>%
  filter(potentially_novel_wen_mcps == "TRUE")
dim(overall_exc_novel)
dim(overall_exc_novel2)
dnmt3a_novel2 <- dnmt3a_novel2 %>%
  filter(potentially_novel_wen_mcps == "TRUE")
dim(dnmt3a_novel)
dim(dnmt3a_novel2)
tet2_novel2 <- tet2_novel2 %>%
  filter(potentially_novel_wen_mcps == "TRUE")
dim(tet2_novel)
dim(tet2_novel2)
asxl1_novel2 <- asxl1_novel2 %>%
  filter(potentially_novel_wen_mcps == "TRUE")
dim(asxl1_novel)
dim(asxl1_novel2)

### Liftover wen_meta from grch38 to grch37:
wen_meta <- wen_meta %>% drop_na()
wen_meta2 <- wen_meta %>%
  rename(CHR = chrom,
         BP = genpos)
wen_meta2 <- unite(wen_meta2, "SNP", CHR:BP, remove = FALSE)
wen_meta3 <- data.table(wen_meta2)
wen_meta3 <- liftover(wen_meta3, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
wen_meta3 <- as_tibble(wen_meta3)

### Separate wen_meta3 by CH type:
wen_meta_dnmt3a <- wen_meta3 %>%
  filter(CH %in% c("All CH", "DNMT3A"))
wen_meta_tet2 <- wen_meta3 %>%
  filter(CH %in% c("All CH", "TET2"))
wen_meta_asxl1 <- wen_meta3 %>%
  filter(CH %in% c("All CH", "ASXL1"))

### Identify potentially novel loci overlapping with wen_meta:
overall_inc_novel3 <- overall_inc_novel2 %>%
  rowwise() %>%
  mutate(min_distance_wen_meta = if_else(any(wen_meta3$CHR == chr),
                                         min(abs(wen_meta3$BP[wen_meta3$CHR == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_meta = if_else(any(wen_meta3$CHR == chr & wen_meta3$BP == pos), "yes", "no"),
         potentially_novel_wen_meta = if_else(!is.na(min_distance_wen_meta) & min_distance_wen_meta > 500000 & overlap_wen_meta == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, min_distance_wen_meta, overlap_wen_meta, potentially_novel_wen_meta, everything())
overall_inc_novel3

overall_exc_novel3 <- overall_exc_novel2 %>%
  rowwise() %>%
  mutate(min_distance_wen_meta = if_else(any(wen_meta3$CHR == chr),
                                         min(abs(wen_meta3$BP[wen_meta3$CHR == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_meta = if_else(any(wen_meta3$CHR == chr & wen_meta3$BP == pos), "yes", "no"),
         potentially_novel_wen_meta = if_else(!is.na(min_distance_wen_meta) & min_distance_wen_meta > 500000 & overlap_wen_meta == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, min_distance_wen_meta, overlap_wen_meta, potentially_novel_wen_meta, everything())
overall_exc_novel3

dnmt3a_novel3 <- dnmt3a_novel2 %>%
  rowwise() %>%
  mutate(min_distance_wen_meta = if_else(any(wen_meta_dnmt3a$CHR == chr),
                                         min(abs(wen_meta_dnmt3a$BP[wen_meta_dnmt3a$CHR == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_meta = if_else(any(wen_meta_dnmt3a$CHR == chr & wen_meta_dnmt3a$BP == pos), "yes", "no"),
         potentially_novel_wen_meta = if_else(!is.na(min_distance_wen_meta) & min_distance_wen_meta > 500000 & overlap_wen_meta == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, min_distance_wen_meta, overlap_wen_meta, potentially_novel_wen_meta, everything())
dnmt3a_novel3

tet2_novel3 <- tet2_novel2 %>%
  rowwise() %>%
  mutate(min_distance_wen_meta = if_else(any(wen_meta_tet2$Chr == chr),
                                         min(abs(wen_meta_tet2$Pos[wen_meta_tet2$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_meta = if_else(any(wen_meta_tet2$Chr == chr & wen_meta_tet2$Pos == pos), "yes", "no"),
         potentially_novel_wen_meta = if_else(!is.na(min_distance_wen_meta) & min_distance_wen_meta > 500000 & overlap_wen_meta == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, min_distance_wen_meta, overlap_wen_meta, potentially_novel_wen_meta, everything())
tet2_novel3

asxl1_novel3 <- asxl1_novel2 %>%
  rowwise() %>%
  mutate(min_distance_wen_meta = if_else(any(wen_meta_asxl1$Chr == chr),
                                         min(abs(wen_meta_asxl1$Pos[wen_meta_asxl1$Chr == chr] - pos), na.rm = TRUE),
                                         Inf),
         overlap_wen_meta = if_else(any(wen_meta_asxl1$Chr == chr & wen_meta_asxl1$Pos == pos), "yes", "no"),
         potentially_novel_wen_meta = if_else(!is.na(min_distance_wen_meta) & min_distance_wen_meta > 500000 & overlap_wen_meta == "no", TRUE, FALSE)) %>%
  ungroup() %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, overlap, min_distance, min_distance_wen_meta, overlap_wen_meta, potentially_novel_wen_meta, everything())
asxl1_novel3

### Filter potentially_novel_wen_meta = TRUE:
overall_inc_novel4 <- overall_inc_novel3 %>%
  filter(potentially_novel_wen_meta == "TRUE")
dim(overall_inc_novel3)
dim(overall_inc_novel4)
overall_exc_novel4 <- overall_exc_novel3 %>%
  filter(potentially_novel_wen_meta == "TRUE")
dim(overall_exc_novel3)
dim(overall_exc_novel4)
dnmt3a_novel4 <- dnmt3a_novel3 %>%
  filter(potentially_novel_wen_meta == "TRUE")
dim(dnmt3a_novel3)
dim(dnmt3a_novel4)
tet2_novel4 <- tet2_novel3 %>%  
  filter(potentially_novel_wen_meta == "TRUE")
dim(tet2_novel3)
dim(tet2_novel4)
asxl1_novel4 <- asxl1_novel3 %>%
  filter(potentially_novel_wen_meta == "TRUE")
dim(asxl1_novel3)
dim(asxl1_novel4)  

### Save:
write_csv(overall_inc_novel4, file = "overall_inc_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv")
write_csv(overall_exc_novel4, file = "overall_exc_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv")
write_csv(dnmt3a_novel4, file = "dnmt3a_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv")
write_csv(tet2_novel4, file = "tet2_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv")
write_csv(asxl1_novel4, file = "asxl1_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv")