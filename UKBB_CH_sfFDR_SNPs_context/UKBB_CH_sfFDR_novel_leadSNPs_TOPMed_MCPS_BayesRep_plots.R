### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(ggpubr)

### Load data:
bick_rbf <- vroom("sffdr_novel_kessler_bick_BayesRep.csv", col_select = c(CH_GWAS, rsID, nearestGene, chr, pos, repBF))
wen_rbf <- vroom("sffdr_novel_kessler_wen_mcps_only_BayesRep.csv", col_select = c(CH_GWAS, rsID, nearestGene, chr, pos, repBF))

### Combine data:
bick_rbf2 <- bick_rbf %>%
  rename(RBF_TOPMed = repBF)
wen_rbf2 <- wen_rbf %>%
  rename(RBF_MCPS = repBF) %>%
  select(rsID, RBF_MCPS)
rbf <- bick_rbf2 %>%
  inner_join(wen_rbf2, by = "rsID")

### Format for plotting:
rbf2 <- rbf %>%
  pivot_longer(cols = starts_with("RBF"),
               names_to = "Study",
               values_to = "RBF") %>%
  mutate(Study = recode(Study,
                        RBF_TOPMed = "TOPMed",
                        RBF_MCPS = "MCPS"))

### Prepare for plotting:
rbf_spaced <- rbf2 %>%
  mutate(CH_GWAS = str_pad(CH_GWAS, width = 15, side = "right"),
         rsID_padded = str_pad(rsID, width = 15, side = "right"),
         nearestGene = str_pad(nearestGene, width = 18, side = "right"),
         label = ifelse(!is.na(rsID),
                   paste0(CH_GWAS, " ", rsID_padded, " ", nearestGene), NA))
custom_order <- c("rs62227732", "rs62550973", "rs4709820", "rs4918404", "rs11668706", 
                  "rs116403919", "rs7170049", "rs9535403", "rs10506957", "rs11596235", 
                  "rs7800079", "rs112540634", "rs1317082", "rs12903325", "rs9267820", 
                  "rs1034332", "rs1474547", "rs17758695", "rs4262954", "rs60626639", 
                  "rs2417347", "rs55903902", "rs35150201", "rs1054684", "rs536897754", 
                  "rs9404952", "rs2237154", "rs1062158", "rs342467", "rs113542380", 
                  "rs11674811", "rs1707334", "rs34016663")
rbf_spaced <- rbf_spaced %>%
  mutate(base_rsID = str_extract(rsID, "rs[0-9]+")) %>%
  mutate(label = factor(label, levels = label[match(custom_order, base_rsID)]))
rbf_spaced <- rbf_spaced %>%
  drop_na()

## Create plot and save
rbf_plot <- ggplot(rbf_spaced, aes(x = label, y = log10(RBF), shape = Study, fill = Study)) +
  geom_point(position = position_dodge(width = 0.6),
             size = 3,
             color = "black", 
             stroke = 0.6,
             na.rm = TRUE) +
  scale_shape_manual(values = c("TOPMed" = 21, "MCPS" = 22)) +
  scale_fill_manual(values = c("TOPMed" = "black", "MCPS" = "white")) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "sfFDR-identified novel lead SNPs",
       y = expression(log[10]*"(RBF)"),
       shape = "Replication Cohort",
       fill = "Replication Cohort") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 18, family = "mono", face = "bold"),   
        axis.text.y = element_text(size = 15, family = "mono", face = "bold"),
        axis.title = element_text(size = 18, family = "mono", face = "bold"),         
        legend.text = element_text(size = 18, family = "mono", face = "bold"),      
        legend.title = element_text(size = 18, family = "mono", face = "bold"),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(color = "black", fill = "white", size = 1),
        strip.text = element_text(family = "mono", size = 10, face = "bold")) +
  ylim(-3, 3)
rbf_plot
ggsave(rbf_plot, file = "rbf_plot_points.pdf")
