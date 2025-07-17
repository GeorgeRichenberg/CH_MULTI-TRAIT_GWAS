### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(stringr)
library(ggforce)

#### Load data:
sffdr_bcx <- vroom("sffdr_all_gws_vs_bcx_traits.csv")
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea, oa, z, p, fp, nearestGene))

### Format sffdr_novel:
sffdr_novel2 <- sffdr_novel %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR \\[|\\]"))
rsid <- sffdr_novel2 %>%
  select(CH_GWAS, rsID) %>%
  rename(rsid = rsID)

#### Filter novel variants from sffdr_bcx:
sffdr_bcx2 <- sffdr_bcx %>%
  rename(rsid = rsID)
sffdr_bcx3 <- sffdr_bcx2 %>%
  inner_join(rsid, by = c("CH_GWAS", "rsid")) %>%
  select(-matches("(^|_)or($|_)|-?l95cl|-?u95cl|^eaf|_eaf"))

### Format sffdr_bcx:
sffdr_bcx4 <- sffdr_bcx3 %>%
  mutate(flip = beta < 0,
         ea = ifelse(flip, oa, ea),
         oa = ifelse(flip, ea, oa),
         beta = ifelse(flip, -beta, beta),
         beta_wbc = ifelse(flip, -beta_wbc, beta_wbc),
         beta_neu = ifelse(flip, -beta_neu, beta_neu),
         beta_lym = ifelse(flip, -beta_lym, beta_lym),
         beta_mon = ifelse(flip, -beta_mon, beta_mon),
         beta_bas = ifelse(flip, -beta_bas, beta_bas),
         beta_eos = ifelse(flip, -beta_eos, beta_eos),
         beta_rbc = ifelse(flip, -beta_rbc, beta_rbc),
         beta_hgb = ifelse(flip, -beta_hgb, beta_hgb),
         beta_hct = ifelse(flip, -beta_hct, beta_hct),
         beta_mch = ifelse(flip, -beta_mch, beta_mch),
         beta_mcv = ifelse(flip, -beta_mcv, beta_mcv),
         beta_mchc = ifelse(flip, -beta_mchc, beta_mchc),
         beta_rdw = ifelse(flip, -beta_rdw, beta_rdw),
         beta_plt = ifelse(flip, -beta_plt, beta_plt),
         beta_mpv = ifelse(flip, -beta_mpv, beta_mpv)) %>%
  select(-flip)

### Calculate id and z-scores:
sffdr_bcx4 <- sffdr_bcx4 %>%
  mutate(id = paste0(
    str_pad(CH_GWAS, width = 13, side = "right"), "  ",
    str_pad(paste0(rsid, "-", ea), width = 14, side = "right"), "  ",
    "chr", str_pad(chr, width = 2, side = "left"), ":", 
    str_pad(pos, width = 9, side = "left"), "  ",
    str_pad(nearestGene, width = 18, side = "right"))) %>%
  rename(pval = p, pval_f = fp)
sffdr_bcx5 <- sffdr_bcx4 %>%
  mutate(z = beta / se,
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
  select(id, CH_GWAS, starts_with("z"), starts_with("pval")) %>%
  select(everything(), -c(z, pval, pval_f))

### Format data for plotting:
trait_order <- c("WBC", "NEU", "LYM", "MON", "BAS", "EOS",
                 "RBC", "HGB", "HCT", "MCH", "MCV", "MCHC", "RDW", 
                 "PLT", "MPV")
sffdr_bcx6 <- sffdr_bcx5 %>%
  pivot_longer(cols = -id,
               names_to = c(".value", "Trait"),
               names_pattern = "^(pval|z)_(.*)$") %>%
  rename(rsID = id, P_value = pval, Z_score = z) %>%
  mutate(Trait = str_to_upper(Trait),
         Trait = factor(Trait, levels = trait_order),
         Significant = P_value < 5e-8,
         Nominally_significant = P_value < 0.05 & P_value >= 5e-8,
         Not_significant = P_value >= 0.05,
         `Z-score direction` = ifelse(Z_score > 0, "Positive", "Negative")) %>%
  drop_na()

### Order y-axis before plotting:
custom_order <- c("rs62227732", "rs62550973", "rs4709820", "rs4918404", "rs11668706", 
                  "rs116403919", "rs7170049", "rs9535403", "rs10506957", "rs11596235", 
                  "rs7800079", "rs112540634", "rs1317082", "rs12903325", "rs9267820", 
                  "rs1034332", "rs1474547", "rs17758695", "rs4262954", "rs60626639", 
                  "rs2417347", "rs55903902", "rs35150201", "rs1054684", "rs536897754", 
                  "rs9404952", "rs2237154", "rs1062158", "rs342467", "rs113542380", 
                  "rs11674811", "rs1707334", "rs34016663")
sffdr_bcx6 <- sffdr_bcx6 %>%
  mutate(base_rsID = str_extract(rsID, "rs[0-9]+"),
         rsID = factor(rsID, levels = rsID[match(custom_order, base_rsID)]))

### Create bubble plot:
plot <- ggplot(sffdr_bcx6, aes(x = Trait, y = rsID, size = Z_score, color = `Z-score direction`)) +
  geom_point(alpha = 0.8) +
  geom_point(data = filter(sffdr_bcx6, Nominally_significant),
             aes(x = Trait, y = rsID, size = Z_score),
             shape = 21, stroke = 1.0, color = "black", fill = NA) +
  geom_point(data = filter(sffdr_bcx6, Significant),
             aes(x = Trait, y = rsID, size = Z_score),
             shape = 21, stroke = 2.0, color = "black", fill = NA) +
  scale_size(range = c(-1, 7.5),
             breaks = seq(-15, 30, by = 15),  
             limits = c(-20, 50),         
             name = "Z-score") +    
  scale_color_manual(values = c("Positive" = "tomato1", "Negative" = "skyblue3"),
                     guide = guide_legend(override.aes = list(size = 8))) +
  theme_minimal() +
  labs(x = "Phenotype",
       y = "sfFDR-identified novel lead SNPs",
       size = "Z-score",
       color = "Z-score") +
  theme(panel.background = element_rect(fill = "white", color = NA),
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
ggsave(plot, file = "ch_bcx_novel_lead_snps_bubble_plot.pdf")