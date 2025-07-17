### Set wd and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(stringr)
library(ggforce)

#### Load data:
sffdr_novel <- vroom("ch_kessler_sffdr_fuma_novel_leadsnps_with_kessler_sumstats.csv", col_select = c(CH_GWAS, rsID, nearestGene))

sffdr_blood <-  vroom("sffdr_novel_finngen_aml_myelodys_myeloprof_myeloid.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea,
                                                                                               beta_aml, se_aml, p_aml,
                                                                                               beta_myelodys, se_myelodys, p_myelodys,
                                                                                               beta_myeloprof, se_myeloprof, p_myeloprof, 
                                                                                               beta_myeloid, se_myeloid, p_myeloid))
sffdr_gene <- vroom("sffdr_novel_plasma_proteins_cebpa_flt3_flt3lg_kit_kitlg_npm1_tp53.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea,
                                                                                                            beta_cebpa, se_cebpa, p_cebpa,
                                                                                                            beta_flt3, se_flt3, p_flt3,
                                                                                                            beta_flt3lg, se_flt3lg, p_flt3lg, 
                                                                                                            beta_kit, se_kit, p_kit,
                                                                                                            beta_kitlg, se_kitlg, p_kitlg,
                                                                                                            beta_npm1, se_npm1, p_npm1,
                                                                                                            beta_tp53, se_tp53, p_tp53))
sffdr_ltl_pc <- vroom("sffdr_novel_ltl.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea,
                                                            beta_ltl_pc1, se_ltl_pc1, pval_ltl_pc1,
                                                            beta_ltl_pc2, se_ltl_pc2, pval_ltl_pc2))
sffdr_mos <- vroom("sffdr_novel_mca_mlox_mloy.csv", col_select = c(CH_GWAS, rsID, chr, pos, ea,
                                                                   beta_mca, se_mca, pval_mca,
                                                                   beta_mlox, se_mlox, pval_mlox,
                                                                   beta_mloy, se_mloy, pval_mloy))

###  Create id column: 
sffdr_blood <- sffdr_blood %>%
  mutate(id = paste(CH_GWAS, rsID, chr, pos, ea, sep = "_")) %>%
  select(-rsID, -chr, -pos, -ea)

sffdr_gene <- sffdr_gene %>%
  mutate(id = paste(CH_GWAS, rsID, chr, pos, ea, sep = "_")) %>%
  select(-rsID, -chr, -pos, -ea)

sffdr_ltl_pc <- sffdr_ltl_pc %>%
  mutate(id = paste(CH_GWAS, rsID, chr, pos, ea, sep = "_")) %>%
  select(-rsID, -chr, -pos, -ea) %>%
  rename_with(.cols = contains("pval_"), .fn = ~ str_replace(., "pval_", "p_"))

sffdr_mos <- sffdr_mos %>%
  mutate(id = paste(CH_GWAS, rsID, chr, pos, ea, sep = "_")) %>%
  select(-rsID, -chr, -pos, -ea) %>%
  rename_with(.cols = contains("pval_"), .fn = ~ str_replace(., "pval_", "p_"))

#### Merge data:
sffdr_all <- sffdr_blood %>%
  inner_join(sffdr_gene, by = "id") %>%
  inner_join(sffdr_ltl_pc, by = "id") %>%
  inner_join(sffdr_mos, by = "id") %>%
  select(id, everything())

### Calculate z-scores:
sffdr_all2 <- sffdr_all %>%
  mutate(across(.cols = matches("^beta_"),
                .fns = ~ .x / get(gsub("^beta_", "se_", cur_column())),
                .names = "z_{.col}")) %>%
  select(id, starts_with("z"), starts_with("p"))

### Format data for plotting:
trait_name_map <- c("aml" = "AML", "myelodys" = "MDS", "myeloprof" = "MPN", "myeloid" = "ML",
                    "cebpa" = "CEBPA", "flt3" = "FLT3", "flt3lg" = "FLT3LG", "kit" = "KIT", "kitlg" = "KITLG", "npm1" = "NPM1", "tp53" = "TP53",
                    "ltl_pc1" = "LTL (PC1)", "ltl_pc2" = "LTL (PC2)",
                    "mca" = "mCA", "mlox" = "mLOX", "mloy" = "mLOY")
trait_order <- unname(trait_name_map)

z_scores <- sffdr_all2 %>%
  select(id, matches("^z_beta_")) %>%
  pivot_longer(cols = -id, names_to = "Trait", values_to = "Z_score") %>%
  mutate(Trait = str_remove(Trait, "^z_beta_"),
         Trait = trait_name_map[Trait])

p_values <- sffdr_all2 %>%
  select(id, matches("^p_")) %>%
  pivot_longer(cols = -id, names_to = "Trait", values_to = "P_value") %>%
  mutate(Trait = str_remove(Trait, "^p_"),
         Trait = trait_name_map[Trait])

sffdr_all3 <- z_scores %>%
  inner_join(p_values, by = c("id", "Trait")) %>%
  rename(rsID = id) %>%
  mutate(Trait = factor(Trait, levels = trait_order),
         Significant = P_value < 5e-8,
         Nominally_significant = P_value < 0.05 & P_value >= 5e-8,
         Not_significant = P_value >= 0.05,
         Direction = ifelse(Z_score > 0, "Positive", "Negative"))
sffdr_all3 <- sffdr_all3 %>%
  mutate(rsID = as.character(rsID))
sffdr_all4 <- sffdr_all3 %>%
  mutate(extracted_rsID = str_extract(rsID, "rs\\d+")) %>%
  left_join(sffdr_novel, by = c("extracted_rsID" = "rsID")) %>%
  mutate(rsID = paste(rsID, nearestGene, sep = "_")) %>%
  select(-extracted_rsID)

### Create id column: 
sffdr_all5 <- sffdr_all4 %>%
  separate(rsID, into = c("ch", "rsid", "chr", "pos", "ea", "gene"), sep = "_")
sffdr_all5 <- sffdr_all5 %>%
  mutate(id = paste0(str_pad(ch, width = 13, side = "right"), "  ",
                     str_pad(paste0(rsid, "-", ea), width = 14, side = "right"), "  ", "chr",
                     str_pad(chr, width = 2, side = "left"), ":", 
                     str_pad(pos, width = 9, side = "left"), "  ",
                     str_pad(gene, width = 18, side = "right"))) %>%
  select(id, everything())

### Order y-axis before plotting:
custom_order <- c("rs62227732", "rs62550973", "rs4709820", "rs4918404", "rs11668706", 
                  "rs116403919", "rs7170049", "rs9535403", "rs10506957", "rs11596235", 
                  "rs7800079", "rs112540634", "rs1317082", "rs12903325", "rs9267820", 
                  "rs1034332", "rs1474547", "rs17758695", "rs4262954", "rs60626639", 
                  "rs2417347", "rs55903902", "rs35150201", "rs1054684", "rs536897754", 
                  "rs9404952", "rs2237154", "rs1062158", "rs342467", "rs113542380", 
                  "rs11674811", "rs1707334", "rs34016663")
sffdr_all5 <- sffdr_all5 %>%
  mutate(base_rsID = str_extract(id, "rs[0-9]+"),
         id = factor(id, levels = id[match(custom_order, base_rsID)]))

### Create bubble plot:
plot <- ggplot(sffdr_all5, aes(x = Trait, y = id, size = Z_score, color = Direction)) +
  geom_point(alpha = 0.8, na.rm = FALSE) +
  geom_point(data = filter(sffdr_all5, Nominally_significant),
             aes(x = Trait, y = id, size = Z_score),
             shape = 21, stroke = 1.0, color = "black", fill = NA, na.rm = FALSE) +
  geom_point(data = filter(sffdr_all5, Significant),
             aes(x = Trait, y = id, size = Z_score),
             shape = 21, stroke = 2.0, color = "black", fill = NA, na.rm = FALSE) +
  scale_size(range = c(-1, 7.5),
             breaks = seq(-15, 30, by = 15),  
             limits = c(-20, 55),         
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
        axis.text.y = element_text(family = "mono", face = "bold", size = 15),
        axis.text.x = element_text(family = "mono", face = "bold", angle = 45, size = 15, vjust = 1.05, hjust = 0.9),
        legend.text = element_text(size = 14, family = "mono", face = "bold"),
        legend.title = element_text(size = 16, family = "mono", face = "bold"))
plot

### Save:
ggsave(plot, file = "ch_bcx_novel_lead_snps_pheno_protein_mca_bubble_plot.pdf")