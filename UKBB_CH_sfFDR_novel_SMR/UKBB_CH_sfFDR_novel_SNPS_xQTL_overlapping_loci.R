### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(readxl)
library(data.table)
library(MungeSumstats)
library(writexl)

### Load data:
ch_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, chr, pos, rsID))
overall_inc_smr <- read_excel("overall_inc_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
overall_exc_smr <- read_excel("overall_exc_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
dnmt3a_smr <- read_excel("dnmt3a_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
tet2_smr <- read_excel("tet2_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
asxl1_smr <- read_excel("asxl1_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")

### Format ch_novel:
ch_novel2 <- ch_novel %>%
  mutate(CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[Overall \\(inc\\)\\]", "Overall (inc)"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[Overall \\(exc\\)\\]", "Overall (exc)"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[DNMT3A\\]", "DNMT3A"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[TET2\\]", "TET2"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[ASXL1\\]", "ASXL1"))
ch_novel2 <- ch_novel2 %>%
  mutate(start = pos - 1000000, 
         end = pos + 1000000)

### Lift smr results from GRCh38 to GRCh37:
overall_inc_smr2 <- overall_inc_smr %>%
  select(-c(chr, start, end, strand, GWAS_LOCUS, Lead_SNP, Lead_SNP_BP)) %>%
  rename(chr = topSNP_chr,
         bp = topSNP_bp)
overall_inc_smr2 <- unite(overall_inc_smr2, "SNP", chr:bp, remove = FALSE)
overall_inc_smr2 <- data.table(overall_inc_smr2)
overall_inc_smr2 <- liftover(overall_inc_smr2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
overall_inc_smr2 <- as_tibble(overall_inc_smr2)
overall_inc_smr2 <- select(overall_inc_smr2, -SNP)
overall_inc_smr2 <- overall_inc_smr2 %>%
  rename(chr = CHR,
         bp = BP)
overall_inc_smr2 <- overall_inc_smr2 %>%
  mutate(chr = as.numeric(chr),
         bp = as.numeric(bp)) %>%
  mutate(CH_GWAS = "Overall (inc)")

overall_exc_smr2 <- overall_exc_smr %>%
  select(-c(chr, start, end, strand, GWAS_LOCUS, Lead_SNP, Lead_SNP_BP)) %>%
  rename(chr = topSNP_chr,
         bp = topSNP_bp)
overall_exc_smr2 <- unite(overall_exc_smr2, "SNP", chr:bp, remove = FALSE)
overall_exc_smr2 <- data.table(overall_exc_smr2)
overall_exc_smr2 <- liftover(overall_exc_smr2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
overall_exc_smr2 <- as_tibble(overall_exc_smr2)
overall_exc_smr2 <- select(overall_exc_smr2, -SNP)
overall_exc_smr2 <- overall_exc_smr2 %>%
  rename(chr = CHR,
         bp = BP)
overall_exc_smr2 <- overall_exc_smr2 %>%
  mutate(chr = as.numeric(chr),
         bp = as.numeric(bp)) %>%
  mutate(CH_GWAS = "Overall (exc)")

dnmt3a_smr2 <- dnmt3a_smr %>%
  select(-c(chr, start, end, strand, GWAS_LOCUS, Lead_SNP, Lead_SNP_BP)) %>%
  rename(chr = topSNP_chr,
         bp = topSNP_bp)
dnmt3a_smr2 <- unite(dnmt3a_smr2, "SNP", chr:bp, remove = FALSE)
dnmt3a_smr2 <- data.table(dnmt3a_smr2)
dnmt3a_smr2 <- liftover(dnmt3a_smr2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
dnmt3a_smr2 <- as_tibble(dnmt3a_smr2)
dnmt3a_smr2 <- select(dnmt3a_smr2, -SNP)
dnmt3a_smr2 <- dnmt3a_smr2 %>%
  rename(chr = CHR,
         bp = BP)
dnmt3a_smr2 <- dnmt3a_smr2 %>%
  mutate(chr = as.numeric(chr),
         bp = as.numeric(bp)) %>%
  mutate(CH_GWAS = "DNMT3A")

tet2_smr2 <- tet2_smr %>%
  select(-c(chr, start, end, strand, GWAS_LOCUS, Lead_SNP, Lead_SNP_BP)) %>%
  rename(chr = topSNP_chr,
         bp = topSNP_bp)
tet2_smr2 <- unite(tet2_smr2, "SNP", chr:bp, remove = FALSE)
tet2_smr2 <- data.table(tet2_smr2)
tet2_smr2 <- liftover(tet2_smr2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
tet2_smr2 <- as_tibble(tet2_smr2)
tet2_smr2 <- select(tet2_smr2, -SNP)
tet2_smr2 <- tet2_smr2 %>%
  rename(chr = CHR,
         bp = BP)
tet2_smr2 <- tet2_smr2 %>%
  mutate(chr = as.numeric(chr),
         bp = as.numeric(bp)) %>%
  mutate(CH_GWAS = "TET2")

asxl1_smr2 <- asxl1_smr %>%
  select(-c(chr, start, end, strand, GWAS_LOCUS, Lead_SNP, Lead_SNP_BP)) %>%
  rename(chr = topSNP_chr,
         bp = topSNP_bp)
asxl1_smr2 <- unite(asxl1_smr2, "SNP", chr:bp, remove = FALSE)
asxl1_smr2 <- data.table(asxl1_smr2)
asxl1_smr2 <- liftover(asxl1_smr2, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
asxl1_smr2 <- as_tibble(asxl1_smr2)
asxl1_smr2 <- select(asxl1_smr2, -SNP)
asxl1_smr2 <- asxl1_smr2 %>%
  rename(chr = CHR,
         bp = BP)
asxl1_smr2 <- asxl1_smr2 %>%
  mutate(chr = as.numeric(chr),
         bp = as.numeric(bp)) %>%
  mutate(CH_GWAS = "ASXL1")

### Split smr results by qtl type:
overall_inc_smr_eqtl <- overall_inc_smr2 %>%
  filter(str_detect(qtl_name, "eQTL"))
overall_inc_smr_pqtl <- overall_inc_smr2 %>%
  filter(str_detect(qtl_name, "pQTL"))
overall_inc_smr_mqtl <- overall_inc_smr2 %>%
  filter(str_detect(qtl_name, "mQTL"))

overall_exc_smr_eqtl <- overall_exc_smr2 %>%
  filter(str_detect(qtl_name, "eQTL"))
overall_exc_smr_pqtl <- overall_exc_smr2 %>%
  filter(str_detect(qtl_name, "pQTL"))
overall_exc_smr_mqtl <- overall_exc_smr2 %>%
  filter(str_detect(qtl_name, "mQTL"))

dnmt3a_smr_eqtl <- dnmt3a_smr2 %>%
  filter(str_detect(qtl_name, "eQTL"))
dnmt3a_smr_pqtl <- dnmt3a_smr2 %>%
  filter(str_detect(qtl_name, "pQTL"))
dnmt3a_smr_mqtl <- dnmt3a_smr2 %>%
  filter(str_detect(qtl_name, "mQTL"))

tet2_smr_eqtl <- tet2_smr2 %>%
  filter(str_detect(qtl_name, "eQTL"))
tet2_smr_pqtl <- tet2_smr2 %>%
  filter(str_detect(qtl_name, "pQTL"))
tet2_smr_mqtl <- tet2_smr2 %>%
  filter(str_detect(qtl_name, "mQTL"))

asxl1_smr_eqtl <- asxl1_smr2 %>%
  filter(str_detect(qtl_name, "eQTL"))
asxl1_smr_pqtl <- asxl1_smr2 %>%
  filter(str_detect(qtl_name, "pQTL"))
asxl1_smr_mqtl <- asxl1_smr2 %>%
  filter(str_detect(qtl_name, "mQTL"))

### Combine data by qtl type:
smr_eqtl <- rbind(overall_inc_smr_eqtl, overall_exc_smr_eqtl, dnmt3a_smr_eqtl, tet2_smr_eqtl, asxl1_smr_eqtl) %>%
  select(CH_GWAS, everything())
dim(smr_eqtl)
smr_pqtl <- rbind(overall_inc_smr_pqtl, overall_exc_smr_pqtl, dnmt3a_smr_pqtl, tet2_smr_pqtl, asxl1_smr_pqtl) %>%
  select(CH_GWAS, everything())
dim(smr_pqtl)
smr_mqtl <- rbind(overall_inc_smr_mqtl, overall_exc_smr_mqtl, dnmt3a_smr_mqtl, tet2_smr_mqtl, asxl1_smr_mqtl) %>%
  select(CH_GWAS, everything())
dim(smr_mqtl)

### Identify smr results overlapping with novel sffdr loci:
ch_novel_eqtl <- ch_novel2 %>%
  inner_join(smr_eqtl, by = c("chr"))
dim(ch_novel_eqtl)
ch_novel_eqtl2 <- ch_novel_eqtl %>%
  group_by(chr) %>%
  filter(bp >= start & bp <= end) %>%
  ungroup()
ch_novel_eqtl2

ch_novel_pqtl <- ch_novel2 %>%
  inner_join(smr_pqtl, by = c("chr"))
dim(ch_novel_pqtl)
ch_novel_pqtl2 <- ch_novel_pqtl %>%
  group_by(chr) %>%
  filter(bp >= start & bp <= end) %>%
  ungroup()
ch_novel_pqtl2

ch_novel_mqtl <- ch_novel2 %>%
  inner_join(smr_mqtl, by = c("chr"))
dim(ch_novel_mqtl)
ch_novel_mqtl2 <- ch_novel_mqtl %>%
  group_by(chr) %>%
  filter(bp >= start & bp <= end) %>%
  ungroup()
ch_novel_mqtl2

### Format:
ch_novel_eqtl2 <- ch_novel_eqtl2 %>%
  select(CH_GWAS.x, rsID, chr, pos, start, end, CH_GWAS.y, qtl_name, probeID, ProbeChr, Probe_bp, topSNP, bp, A1, A2, Freq, b_GWAS, se_GWAS, p_GWAS, b_eQTL, se_eQTL, p_eQTL, b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI, gene_id, index)
ch_novel_pqtl2 <- ch_novel_pqtl2 %>%
  select(CH_GWAS.x, rsID, chr, pos, start, end, CH_GWAS.y, qtl_name, probeID, ProbeChr, Probe_bp, topSNP, bp, A1, A2, Freq, b_GWAS, se_GWAS, p_GWAS, b_eQTL, se_eQTL, p_eQTL, b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI, gene_id, index)
ch_novel_mqtl2 <- ch_novel_mqtl2 %>%
  select(CH_GWAS.x, rsID, chr, pos, start, end, CH_GWAS.y, qtl_name, probeID, ProbeChr, Probe_bp, topSNP, bp, A1, A2, Freq, b_GWAS, se_GWAS, p_GWAS, b_eQTL, se_eQTL, p_eQTL, b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI, gene_id, index)

### Save:
write_csv(ch_novel_eqtl2, file = "crosstrait_ch_novel_eqtl_p<1e-3_heidi>0.05.csv")
write_csv(ch_novel_pqtl2, file = "crosstrait_ch_novel_pqtl_p<1e-3_heidi>0.05.csv")
write_csv(ch_novel_mqtl2, file = "crosstrait_ch_novel_mqtl_p<1e-3_heidi>0.05.csv")