### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv")
dnmt3a <- vroom("dnmt3a_kessler_sffdr_fp.txt.gz")
tet2 <- vroom("tet2_kessler_sffdr_fp.txt.gz")
asxl1 <- vroom("asxl1_kessler_sffdr_fp.txt.gz")

### Filter novel overall leadsnps:
overall_inc_novel <- novel %>%
  filter(CH_GWAS == "sfFDR [Overall (inc)]")
overall_exc_novel <- novel %>%
  filter(CH_GWAS == "sfFDR [Overall (exc)]")

### Create uniqid in ch subtype data:
dnmt3a <- dnmt3a %>%
  mutate(alleles = pmap_chr(list(A1, A2), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)
tet2 <- tet2 %>%
  mutate(alleles = pmap_chr(list(A1, A2), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)
asxl1 <- asxl1 %>%
  mutate(alleles = pmap_chr(list(A1, A2), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)

### Rename p-vlaue column:
dnmt3a <- dnmt3a %>%
  rename(chr = CHR,
         pos = BP,
         ea = A1,
         oa = A2,
         eaf = EAF,
         beta = Beta,
         se = SE,
         p = P,
         fp = FP)
tet2 <- tet2 %>%
  rename(chr = CHR,
         pos = BP,
         ea = A1,
         oa = A2,
         eaf = EAF,
         beta = Beta,
         se = SE,
         p = P,
         fp = FP)
asxl1 <- asxl1 %>%
  rename(chr = CHR,
         pos = BP,
         ea = A1,
         oa = A2,
         eaf = EAF,
         beta = Beta,
         se = SE,
         p = P,
         fp = FP)

### Merge CH datasets with novel overall_inc SNPs:
overall_inc_novel2 <- overall_inc_novel %>%
  select(uniqID)
dnmt3a_inc <- overall_inc_novel2 %>%
  inner_join(dnmt3a, by = "uniqID")
dim(overall_inc_novel2)
dim(dnmt3a_inc)
tet2_inc <- overall_inc_novel2 %>%
  inner_join(tet2, by = "uniqID")
dim(overall_inc_novel2)
dim(tet2_inc)
asxl1_inc <- overall_inc_novel2 %>%
  inner_join(asxl1, by = "uniqID")
dim(overall_inc_novel2)
dim(asxl1_inc)

### Merge CH datasets with novel overall_exc SNPs:
overall_exc_novel2 <- overall_exc_novel %>%
  select(uniqID)
dnmt3a_exc <- overall_exc_novel2 %>%
  inner_join(dnmt3a, by = "uniqID")
dim(overall_exc_novel2)
dim(dnmt3a_exc)
tet2_exc <- overall_exc_novel2 %>%
  inner_join(tet2, by = "uniqID")
dim(overall_exc_novel2)
dim(tet2_exc)
asxl1_exc <- overall_exc_novel2 %>%
  inner_join(asxl1, by = "uniqID")
dim(overall_exc_novel2)
dim(asxl1_exc)

### Provide trait names per dataset:
dnmt3a_inc <- dnmt3a_inc %>%
  mutate(CH_GWAS = "Overall (inc)",
         CH_subtype_GWAS = "DNMT3A") %>%
  select(CH_GWAS, CH_subtype_GWAS, SNP, chr, pos, uniqID, ea, oa, eaf, beta, se, p, fp)
tet2_inc <- tet2_inc %>%
  mutate(CH_GWAS = "Overall (inc)",
         CH_subtype_GWAS = "TET2") %>%
  select(CH_GWAS, CH_subtype_GWAS, SNP, chr, pos, uniqID, ea, oa, eaf, beta, se, p, fp)
asxl1_inc <- asxl1_inc %>%
  mutate(CH_GWAS = "Overall (inc)",
         CH_subtype_GWAS = "ASXL1") %>%
  select(CH_GWAS, CH_subtype_GWAS, SNP, chr, pos, uniqID, ea, oa, eaf, beta, se, p, fp)

dnmt3a_exc <- dnmt3a_exc %>%
  mutate(CH_GWAS = "Overall (exc)",
         CH_subtype_GWAS = "DNMT3A") %>%
  select(CH_GWAS, CH_subtype_GWAS, SNP, chr, pos, uniqID, ea, oa, eaf, beta, se, p, fp)
tet2_exc <- tet2_exc %>%
  mutate(CH_GWAS = "Overall (exc)",
         CH_subtype_GWAS = "TET2") %>%
  select(CH_GWAS, CH_subtype_GWAS, SNP, chr, pos, uniqID, ea, oa, eaf, beta, se, p, fp)
asxl1_exc <- asxl1_exc %>%
  mutate(CH_GWAS = "Overall (exc)",
         CH_subtype_GWAS = "ASXL1") %>%
  select(CH_GWAS, CH_subtype_GWAS, SNP, chr, pos, uniqID, ea, oa, eaf, beta, se, p, fp)

### Merge datasets and combine with overall_inc_novel + overall_exc_novel:
ch_overall <- rbind(dnmt3a_inc, tet2_inc, asxl1_inc,
                    dnmt3a_exc, tet2_exc, asxl1_exc)
ch_overall <- ch_overall %>%
  rename(rsID = SNP) %>%
  mutate(z = beta / se,
         or = exp(beta),  
         "95_lci" = exp(beta - 1.96 * se),
         "95_uci" = exp(beta + 1.96 * se))
overall_inc_novel3 <- overall_inc_novel %>%
  select(-c(GenomicLocus, nIndSigSNPs))

desired_order <- c("A", "DNMT3A", "TET2", "ASXL1")
ch_overall2 <- bind_rows(ch_overall, overall_inc_novel3, overall_exc_novel3)
ch_overall3 <- ch_overall2 %>%
  select(CH_GWAS, CH_subtype_GWAS, rsID, chr, pos, uniqID, ea, oa, eaf, beta, se, z, or, `95_lci`, `95_uci`, p, fp, nearestGene) %>%
  mutate(CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[Overall \\(inc\\)\\]", "Overall (inc)"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[Overall \\(exc\\)\\]", "Overall (exc)")) %>%
  mutate(CH_subtype_GWAS = replace_na(CH_subtype_GWAS, "A")) %>%
  mutate(CH_subtype_GWAS = factor(CH_subtype_GWAS, levels = desired_order, ordered = TRUE)) %>%
  arrange(chr, pos, CH_subtype_GWAS) %>%
  mutate(CH_subtype_GWAS = as.character(CH_subtype_GWAS))

### Save:
write_csv(ch_overall3, file = "novel_overall_ch_leadsnps_vs_subtype_specific_ch.csv")