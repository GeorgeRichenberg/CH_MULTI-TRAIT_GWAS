### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)

### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv")
overall_afr <- vroom("overall_inc_kessler_afr_GCST90165263_novel_snps.csv", col_select = -c(chr_afr, bp_afr))
overall_eas <- vroom("overall_inc_kessler_eas_GCST90165265_novel_snps.csv", col_select = -c(chr_eas, bp_eas))
overall_sas <- vroom("overall_inc_kessler_sas_GCST90165269_novel_snps.csv", col_select = -c(chr_eas, bp_eas))

### Format afR, eas, sas:
overall_afr <- overall_afr %>%
  rename(rsID = rsid) %>%
  mutate(or_afr = exp(beta_afr),
         z_afr = beta_afr / se_afr,
         l95cl_afr = exp(beta_afr - 1.96 * se_afr),
         u95cl_afr = exp(beta_afr + 1.96 * se_afr)) %>%
  select(rsID, ea_afr, oa_afr, eaf_afr, beta_afr, se_afr, z_afr, or_afr, `l95cl_afr`, `u95cl_afr`, p_afr)
overall_eas <- overall_eas %>%
  rename(rsID = rsid) %>%
  mutate(or_eas = exp(beta_eas),
         z_eas = beta_eas / se_eas,
         l95cl_eas = exp(beta_eas - 1.96 * se_eas),
         u95cl_eas = exp(beta_eas + 1.96 * se_eas)) %>%
  select(rsID, ea_eas, oa_eas, eaf_eas, beta_eas, se_eas, z_eas, or_eas, `l95cl_eas`, `u95cl_eas`, p_eas)
overall_sas <- overall_sas %>%
  rename(rsID = rsid,
         ea_sas = ea_eas,
         oa_sas = oa_eas,
         eaf_sas = eaf_eas,
         beta_sas = beta_eas, 
         se_sas = se_eas,
         p_sas = p_eas) %>%
  mutate(or_sas = exp(beta_sas),
         z_sas = beta_sas / se_sas,
         l95cl_sas = exp(beta_sas - 1.96 * se_sas),
         u95cl_sas = exp(beta_sas + 1.96 * se_sas)) %>%
  select(rsID, ea_sas, oa_sas, eaf_sas, beta_sas, se_sas, z_sas, or_sas, `l95cl_sas`, `u95cl_sas`, p_sas)

### Format sffdr novel:
sffdr_novel <- sffdr_novel %>%
  mutate(CH_GWAS = CH_GWAS %>%
           str_remove_all("sfFDR\\s*") %>%   
           str_remove_all("\\[|\\]")) %>%
  select(-c(GenomicLocus, nIndSigSNPs, uniqID))

### Merge:
sffdr_novel2 <- sffdr_novel %>%
  left_join(overall_afr, by = "rsID") %>%
  left_join(overall_eas, by = "rsID") %>%
  left_join(overall_sas, by = "rsID")

### Save:
write_csv(sffdr_novel2, file = "sffdr_novel_afr_eas_sas.csv")

### Format data for plotting: 
gene_order <- c("ZDHHC18", "PIK3R3", "HEATR5B", "THADA", "AFF1", "NDFIP1", "JARID2",
                "HLA-G", "TMPOP1", "HLA-DOB", "AC093107.7", "MYRF", "ATF7IP", "CTCF",
                "RP4-536B24.3", "BCL2", "NRIP1:AF127577.11", "LINC00649")
sffdr_novel3 <- sffdr_novel2 %>%
  filter(CH_GWAS %in% "Overall (inc)") %>%
  mutate(id = paste0(str_pad(CH_GWAS, width = 13, side = "right"), "  ",
                     str_pad(paste0(rsID, "-", ea), width = 13, side = "right"), "  ", "chr", 
                     str_pad(chr, width = 2, side = "left"), ":", 
                     str_pad(pos, width = 9, side = "left"), "  ",
                     str_pad(nearestGene, width = 17, side = "right")))
sffdr_long <- sffdr_novel3 %>%
  rename(l95cl = `95_lci`, 
         u95cl = `95_uci`) %>%
  select(id, nearestGene, or, l95cl, u95cl,
         or_afr, l95cl_afr, u95cl_afr,
         or_eas, l95cl_eas, u95cl_eas,
         or_sas, l95cl_sas, u95cl_sas) %>%
  pivot_longer(cols = c(or, l95cl, u95cl,
                        or_afr, l95cl_afr, u95cl_afr,
                        or_eas, l95cl_eas, u95cl_eas,
                        or_sas, l95cl_sas, u95cl_sas),
    names_to = c(".value", "Ancestry"),
    names_pattern = "^(or|l95cl|u95cl)_?(.*)$") %>%
  mutate(Ancestry = case_when(Ancestry == "" ~ "EUR",
                              Ancestry == "afr"  ~ "AFR",
                              Ancestry == "eas"  ~ "EAS",
                              Ancestry == "sas"  ~ "SAS", TRUE ~ Ancestry),
    Ancestry = factor(Ancestry, levels = c("SAS", "EAS", "AFR", "EUR")),
    nearestGene = factor(nearestGene, levels = gene_order),
    logOR = log(or),
    logL = log(l95cl),
    logU = log(u95cl)) %>%
  arrange(nearestGene, Ancestry) %>%
  mutate(id = factor(id, levels = rev(unique(id))))

### Create forest plot:
forest_plot_log <- ggplot(sffdr_long, 
  aes(x = id, y = or, ymin = l95cl, ymax = u95cl, color = Ancestry)) +
  geom_pointrange(position = position_dodge(width = 0.6),
                  size = 0.5,
                  shape = 16, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  coord_flip(ylim = c(0, 2)) +
  labs(x = NULL,
       y = expression("OR ("%+-%"95% CI)"),
      color = "Ancestry") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(family = "mono", size = 16, face = "bold"),
        axis.text.y = element_text(family = "mono", size = 14, face = "bold"),
        axis.title.x = element_text(family = "mono", size = 16, face = "bold"),
        legend.text = element_text(family = "mono", size = 16),
        legend.title = element_text(family = "mono", size = 13, face = "bold"))
forest_plot_log

### Save:
ggsave(forest_plot_log, file = "novel_overll_inc_ch_loci_afr_eas_sas_forest_plot.pdf")
