### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load all loci:
all_loci <- vroom("ch_kessler_sffdr_fuma_leadsnps_all_loci.csv")

### Load novel loci:
overall_inc_novel <- vroom("overall_inc_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene))
overall_exc_novel <- vroom("overall_exc_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene))
dnmt3a_novel <- vroom("dnmt3a_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene))
tet2_novel <- vroom("tet2_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene))
asxl1_novel <- vroom("asxl1_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, p, nearestGene))

### Load Kessler data:
overall_inc_kessler <- vroom("overall_GCST90165267_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall))
overall_exc_kessler <- vroom("overall_GCST90165261_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(CHR, BP, ea_overall, oa_overall, eaf_overall, beta_overall, se_overall, pval_overall))
dnmt3a_kessler<- vroom("dnmt3a_GCST90165271_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(CHR, BP, ea_dnmt3a, oa_dnmt3a, eaf_dnmt3a, beta_dnmt3a, se_dnmt3a, pval_dnmt3a))
tet2_kessler <- vroom("tet2_GCST90165281_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(CHR, BP, ea_tet2, oa_tet2, eaf_tet2, beta_tet2, se_tet2, pval_tet2))
asxl1_kessler <- vroom("asxl1_GCST90165259_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz", col_select = c(CHR, BP, ea_asxl1, oa_asxl1, eaf_asxl1, beta_asxl1, se_asxl1, pval_asxl1))

### Rename data:
all_loci <- all_loci %>%
  rename(fp = p)
overall_inc_kessler2 <- overall_inc_kessler %>%
  rename(ea = ea_overall,
         oa = oa_overall,
         eaf = eaf_overall,
         beta = beta_overall,
         se = se_overall,
         p =  pval_overall)
overall_exc_kessler2 <- overall_exc_kessler %>%
  rename(ea = ea_overall,
         oa = oa_overall,
         eaf = eaf_overall,
         beta = beta_overall,
         se = se_overall,
         p =  pval_overall)
dnmt3a_kessler2 <- dnmt3a_kessler %>%
  rename(ea = ea_dnmt3a,
         oa = oa_dnmt3a,
         eaf = eaf_dnmt3a,
         beta = beta_dnmt3a,
         se = se_dnmt3a,
         p =  pval_dnmt3a)
tet2_kessler2 <- tet2_kessler %>%
  rename(ea = ea_tet2,
         oa = oa_tet2,
         eaf = eaf_tet2,
         beta = beta_tet2,
         se = se_tet2,
         p =  pval_tet2)
asxl1_kessler2 <- asxl1_kessler %>%
  rename(ea = ea_asxl1,
         oa = oa_asxl1,
         eaf = eaf_asxl1,
         beta = beta_asxl1,
         se = se_asxl1,
         p =  pval_asxl1)

### Create uniqid in kessler:
overall_inc_kessler2 <- overall_inc_kessler2 %>%
  mutate(alleles = pmap_chr(list(ea, oa), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)
overall_exc_kessler2 <- overall_exc_kessler2 %>%
  mutate(alleles = pmap_chr(list(ea, oa), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)
dnmt3a_kessler2 <- dnmt3a_kessler2 %>%
  mutate(alleles = pmap_chr(list(ea, oa), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)
tet2_kessler2 <- tet2_kessler2 %>%
  mutate(alleles = pmap_chr(list(ea, oa), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)
asxl1_kessler2 <- asxl1_kessler2 %>%
  mutate(alleles = pmap_chr(list(ea, oa), ~paste(sort(c(..1, ..2)), collapse = ":")),
         uniqID = paste(CHR, BP, alleles, sep = ":")) %>%
  select(-alleles)

### Separate all_loci by ch type:
overall_inc_loci <- all_loci %>%
  filter(grepl("Overall \\(inc\\)", CH_GWAS))
overall_inc_loci
overall_exc_loci <- all_loci %>%
  filter(grepl("Overall \\(exc\\)", CH_GWAS))
overall_exc_loci 
dnmt3a_loci <- all_loci %>%
  filter(grepl("DNMT3A", CH_GWAS))
dnmt3a_loci
tet2_loci <- all_loci %>%
  filter(grepl("TET2", CH_GWAS))
tet2_loci
asxl1_loci <- all_loci %>%
  filter(grepl("ASXL1", CH_GWAS))
asxl1_loci

### Join loci with kessler:
overall_inc_loci2 <- overall_inc_loci %>%
  inner_join(overall_inc_kessler2, by = "uniqID") %>%
  select(-c(CHR, BP))
dim(overall_inc_loci)
dim(overall_inc_loci2)
overall_exc_loci2 <- overall_exc_loci %>%
  inner_join(overall_exc_kessler2, by = "uniqID") %>%
  select(-c(CHR, BP))
dim(overall_exc_loci)
dim(overall_exc_loci2)
dnmt3a_loci2 <- dnmt3a_loci %>%
  inner_join(dnmt3a_kessler2, by = "uniqID") %>%
  select(-c(CHR, BP))
dim(dnmt3a_loci)
dim(dnmt3a_loci2)
tet2_loci2 <- tet2_loci %>%
  inner_join(tet2_kessler2, by = "uniqID") %>%
  select(-c(CHR, BP))
dim(tet2_loci)
dim(tet2_loci2)
asxl1_loci2 <- asxl1_loci %>%
  inner_join(asxl1_kessler2, by = "uniqID") %>%
  select(-c(CHR, BP))
dim(asxl1_loci)
dim(asxl1_loci2)

### Combine loci:
all_loci3 <- rbind(overall_inc_loci2,
                   overall_exc_loci2,
                   dnmt3a_loci2,
                   tet2_loci2,
                   asxl1_loci2)
all_loci3 <- all_loci3 %>%
  arrange(chr, pos)

### Calculate z-scores and OR:
all_loci4 <- all_loci3 %>%
  mutate(z = beta / se,
         or = exp(beta),  
         "95_lci" = exp(beta - 1.96 * se),
         "95_uci" = exp(beta + 1.96 * se)) %>%
  select(CH_GWAS, GenomicLocus, nIndSigSNPs, rsID, chr, pos, uniqID, ea, oa, eaf, beta, se, z, or, "95_lci", "95_uci", p, fp, nearestGene)

### Select sfFDR loci only:
all_loci5 <- all_loci4 %>%
  filter(grepl("sfFDR", CH_GWAS))

### Filter novel loci:
novel_loci <- rbind(overall_inc_novel, overall_exc_novel, dnmt3a_novel, tet2_novel, asxl1_novel)
all_loci_novel <- novel_loci %>%
  semi_join(all_loci5, by = c("CH_GWAS", "rsID"))                   

### Save:
write_csv(all_loci4, file = "ch_kessler_sffdr_fuma_all_leadsnps_with_kessler_sumstats.csv")
write_csv(all_loci5, file = "ch_kessler_sffdr_fuma_only_leadsnps_with_kessler_sumstats.csv")
write_csv(all_loci_novel, file = "ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv")
