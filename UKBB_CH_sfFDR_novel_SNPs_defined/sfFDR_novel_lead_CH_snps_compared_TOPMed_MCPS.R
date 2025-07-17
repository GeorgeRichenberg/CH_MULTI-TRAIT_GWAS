### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv")
bick <- vroom("bick_topmed_ch_gwas_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz")
wen_overall <- vroom("wen_overall_GCST90435341_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
wen_dnmt3a <- vroom("wen_dnmt3a_GCST90435342_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
wen_tet2 <- vroom("wen_tet2_GCST90435343_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
wen_asxl1 <- vroom("wen_asxl1_GCST90435344_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Format kessler:
sffdr_novel <- sffdr_novel %>%
  rename(ea_kessler = ea,
         oa_kessler = oa, 
         eaf_kessler = eaf,
         beta_kessler = beta,
         se_kessler = se,
         z_kessler = z,
         or_kessler = or,
         `95_lci_kessler` = `95_lci`,
         `95_uci_kessler` = `95_uci`,
         p_kessler = p,
         fp_sffdr = fp) %>%
  mutate(SNP = str_extract(uniqID, "^[^:]+:[^:]+")) %>% 
  mutate(SNP = str_replace(SNP, ":", "_")) 

### Combine sffdr novel and bick:
novel_sffdr_bick <- sffdr_novel %>%
  inner_join(bick, by = "SNP")

### Calcualte z-scores:
novel_sffdr_bick <- novel_sffdr_bick %>%
  mutate(z_bick = beta_bick / se_bick,
         or_bick = exp(beta_bick),
         `95_lci_bick` = exp(beta_bick - 1.96 * se_bick),
         `95_uci_bick` = exp(beta_bick + 1.96 * se_bick))
novel_sffdr_bick <- novel_sffdr_bick %>%
  select(CH_GWAS, GenomicLocus, chr, pos, rsID,
         ea_kessler, oa_kessler, eaf_kessler, beta_kessler, se_kessler, or_kessler, "95_lci_kessler", "95_uci_kessler", z_kessler, p_kessler, fp_sffdr, nearestGene,
         ea_bick, oa_bick, or_bick, eaf_bick, beta_bick, se_bick, or_bick, "95_lci_bick", "95_uci_bick", z_bick, pval_bick)

### Save:
write_csv(novel_sffdr_bick, file = "sffdr_novel_leadsnps_bick_compared.csv")

### Prepare Wen data:
wen_overall <- wen_overall %>%
  rename(rsID = rsid)
wen_dnmt3a <- wen_dnmt3a %>%
  rename(rsID = rsid)
wen_tet2 <- wen_tet2 %>%
  rename(rsID = rsid)
wen_asxl1 <- wen_asxl1 %>%
  rename(rsID = rsid)

### Combine sffdr novel and wen:
sffdr_novel_overall <- sffdr_novel %>%
  filter(grepl("Overall", CH_GWAS, ignore.case = TRUE))
sffdr_novel_dnmt3a <- sffdr_novel %>%
  filter(grepl("DNMT3A", CH_GWAS, ignore.case = TRUE))
sffdr_novel_tet2 <- sffdr_novel %>%
  filter(grepl("TET2", CH_GWAS, ignore.case = TRUE))
sffdr_novel_asxl1 <- sffdr_novel %>%
  filter(grepl("ASXL1", CH_GWAS, ignore.case = TRUE))

sffdr_novel_overall2 <- sffdr_novel_overall %>%
  inner_join(wen_overall, by = "rsID") %>%
  rename_with(~ gsub("_overall", "", .x, ignore.case = TRUE))
sffdr_novel_overall_dnmt3a <- sffdr_novel_overall %>%
  inner_join(wen_dnmt3a, by = "rsID") %>%
  rename_with(~ gsub("_dnmt3a", "", .x, ignore.case = TRUE)) %>%
  mutate(CH_subtype_GWAS = "DNMT3A")
sffdr_novel_overall_tet2 <- sffdr_novel_overall %>%
  inner_join(wen_tet2, by = "rsID") %>%
  rename_with(~ gsub("_tet2", "", .x, ignore.case = TRUE)) %>%
  mutate(CH_subtype_GWAS = "TET2")
sffdr_novel_overall_asxl1 <- sffdr_novel_overall %>%
  inner_join(wen_asxl1, by = "rsID") %>%
  rename_with(~ gsub("_asxl1", "", .x, ignore.case = TRUE)) %>%
  mutate(CH_subtype_GWAS = "ASXL1")
sffdr_novel_dnmt3a2 <- sffdr_novel_dnmt3a %>%
  inner_join(wen_dnmt3a, by = "rsID") %>%
  rename_with(~ gsub("_dnmt3a", "", .x, ignore.case = TRUE))
sffdr_novel_tet22 <- sffdr_novel_tet2 %>%
  inner_join(wen_tet2, by = "rsID") %>%
  rename_with(~ gsub("_tet2", "", .x, ignore.case = TRUE))
sffdr_novel_asxl12 <- sffdr_novel_asxl1 %>%
  inner_join(wen_asxl1, by = "rsID") %>%
  rename_with(~ gsub("_asxl1", "", .x, ignore.case = TRUE)) 

sffdr_novel_mcps <- bind_rows(sffdr_novel_overall2, sffdr_novel_overall_dnmt3a, sffdr_novel_overall_tet2, sffdr_novel_overall_asxl1, sffdr_novel_dnmt3a2, sffdr_novel_tet22, sffdr_novel_asxl12)

### Calcualte z-scores, or and format:
sffdr_novel_mcps <- sffdr_novel_mcps %>%
  mutate(z_wen = beta_wen / se_wen,
         or_wen = exp(beta_wen),
         `95_lci_wen` = exp(beta_wen - 1.96 * se_wen),
         `95_uci_wen` = exp(beta_wen + 1.96 * se_wen))
sffdr_novel_mcps <- sffdr_novel_mcps %>%
  select(CH_GWAS, CH_subtype_GWAS, GenomicLocus, chr, pos, rsID,
         ea_kessler, oa_kessler, eaf_kessler, beta_kessler, se_kessler, or_kessler, "95_lci_kessler", "95_uci_kessler", z_kessler, p_kessler,
         ea_wen, oa_wen, eaf_wen, beta_wen, se_wen, or_wen, "95_lci_wen", "95_uci_wen", z_wen, pval_wen,
         fp_sffdr, nearestGene)
sffdr_novel_mcps <- sffdr_novel_mcps %>%
  mutate(CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[Overall \\(inc\\)\\]", "Overall (inc)"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[Overall \\(exc\\)\\]", "Overall (exc)"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[DNMT3A\\]", "DNMT3A"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[TET2\\]", "TET2"),
         CH_GWAS = str_replace(CH_GWAS, "sfFDR \\[ASXL1\\]", "ASXL1"))
desired_order <- c("A", "DNMT3A", "TET2", "ASXL1")
sffdr_novel_mcps2 <- sffdr_novel_mcps %>%
  mutate(CH_subtype_GWAS = replace_na(CH_subtype_GWAS, "A")) %>%
  mutate(CH_subtype_GWAS = factor(CH_subtype_GWAS, levels = desired_order, ordered = TRUE)) %>%
  arrange(chr, pos, CH_subtype_GWAS)

### Save:
write_csv(sffdr_novel_mcps2, file = "sffdr_novel_leadsnps_wen_mcps_compared.csv")