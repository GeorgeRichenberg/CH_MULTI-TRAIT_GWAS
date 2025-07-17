### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(stringr)
library(ggforce)

#### Load data:
kessler_lead <- vroom("kessler_lead_snps_15_bcx_traits.csv")

### Calculate z-scores:
kessler_lead <- kessler_lead %>%
  select(rsid, CHR, BP, EA, OA, Beta_UKB, SE_UKB, Pval_UKB, nearestGene,
         beta_wbc, se_wbc, pval_wbc,
         beta_neu, se_neu, pval_neu,
         beta_lym, se_lym, pval_lym,
         beta_mon, se_mon, pval_mon,
         beta_bas, se_bas, pval_bas,
         beta_eos, se_eos, pval_eos,
         beta_rbc, se_rbc, pval_rbc,
         beta_hgb, se_hgb, pval_hgb,
         beta_hct, se_hct, pval_hct,
         beta_mch, se_mch, pval_mch,
         beta_mcv, se_mcv, pval_mcv,
         beta_mchc, se_mchc, pval_mchc,
         beta_rdw, se_rdw, pval_rdw,
         beta_plt, se_plt, pval_plt,
         beta_mpv, se_mpv, pval_mpv)
kessler_lead <- kessler_lead %>%
  mutate(id = paste(rsid, CHR, BP, EA, nearestGene, sep = "_"))
kessler_lead2 <- kessler_lead %>%
  mutate(z_UKB = Beta_UKB / SE_UKB,
         z_wbc = beta_wbc / se_wbc,
         z_neu = beta_neu / se_neu,
         z_lym = beta_lym / se_lym,
         z_mon = beta_mon / se_mon,
         z_bas = beta_bas / se_bas,
         z_eos = beta_bas / se_bas,
         z_rbc = beta_rbc / se_rbc,
         z_hgb = beta_hgb / se_hgb,
         z_hct = beta_hct / se_hct,
         z_mch = beta_mch / se_mch,
         z_mcv = beta_mch / se_mch,
         z_mchc = beta_mchc / se_mchc,
         z_rdw = beta_rdw / se_rdw,
         z_plt = beta_plt / se_plt,
         z_mpv = beta_mpv / se_mpv) %>%
  select(id, starts_with("z"), starts_with("p"))

### Format data for plotting:
trait_order <- c("WBC", "NEU", "LYM", "MON", "BAS", "EOS",
                 "RBC", "HGB", "HCT", "MCH", "MCV", "MCHC", "RDW", 
                 "PLT", "MPV")
kessler_lead3 <- kessler_lead2 %>%
  pivot_longer(
    cols = -id,
    names_to = c(".value", "Trait"),
    names_pattern = "^(pval|z)_(.*)$") %>%
  rename(rsID = id, P_value = pval, Z_score = z) %>%
  mutate(
    Trait = str_to_upper(Trait),
    Trait = factor(Trait, levels = trait_order),
    Significant = P_value < 5e-8,
    Nominally_significant = P_value < 0.05 & P_value >= 5e-8,
    Not_significant = P_value >= 0.05,
    `Z-score direction` = ifelse(Z_score > 0, "Positive", "Negative")) %>%
  drop_na()

### Create id column: 
kessler_lead4 <- kessler_lead3 %>%
  separate(rsID, into = c("rsid", "chr", "pos", "ea", "gene"), sep = "_")
kessler_lead4 <- kessler_lead4 %>%
  mutate(id = paste0("Overall (inc)", "  ", 
         str_pad(paste0(rsid, "-", ea), width = 14, side = "right"), "  ",
         "chr", str_pad(chr, width = 2, side = "left"), ":", 
         str_pad(pos, width = 9, side = "left"), "  ",
         str_pad(gene, width = 18, side = "right"))) %>%
  select(id, everything())

### Order y-axis before plotting:
custom_order <- c(
  "rs361725", "rs62235635", "rs2834708", "rs7272062", "rs8088824",
  "rs80093687", "rs7406485", "rs72698720", "rs16930705", "rs10890839",
  "rs7817503", "rs11762008", "rs2275652", "rs2396004", "rs9264393",
  "rs7705526", "rs144317085", "rs2305407", "rs6440667", "rs62175229",
  "rs3219104")
kessler_lead4 <- kessler_lead4 %>%
  mutate(base_rsID = str_extract(id, "rs[0-9]+"),
         id = factor(id, levels = id[match(custom_order, base_rsID)]))

### Create bubble plot:
plot <- ggplot(kessler_lead4, aes(x = Trait, y = id, size = Z_score, color = `Z-score direction`)) +
  geom_point(alpha = 0.8) +
  geom_point(data = filter(kessler_lead4, Nominally_significant),
             aes(x = Trait, y = id, size = Z_score),
             shape = 21, stroke = 1.0, color = "black", fill = NA) +
  geom_point(data = filter(kessler_lead4, Significant),
             aes(x = Trait, y = id, size = Z_score),
             shape = 21, stroke = 2.0, color = "black", fill = NA) +
  scale_size(range = c(-10, 15),
             breaks = seq(-40, 30, by = 20),  
             limits = c(-100, 100),         
             name = "Z-score") +    
  scale_color_manual(values = c("Positive" = "tomato1", "Negative" = "skyblue3"),
                     guide = guide_legend(override.aes = list(size = 8))) +
  theme_minimal() +
  labs(x = "Phenotype",
       y = "Lead SNPs from Kessler et al. UKBB CH GWAS",
       size = "Z-score",
       color = "Z-score") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.justification = c(0, 0),
    legend.box.just = "left",
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -175),
    axis.title = element_text(size = 18, family = "mono", face = "bold"),
    text = element_text(family = "mono", face = "bold", size = 16),
    legend.text = element_text(size = 14, family = "mono", face = "bold"),
    legend.title = element_text(size = 16, family = "mono", face = "bold"))
plot

### Save:
ggsave(plot, file = "Kessler_ukbb_lead_snps_15_bcx_bubble_plot.pdf")
