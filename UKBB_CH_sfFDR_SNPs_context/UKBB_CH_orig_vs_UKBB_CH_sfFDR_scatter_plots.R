### Set working directory and load required packages:
setwd("/Users/cs21140/Downloads/ch_cross_trait")
library(vroom)
library(tidyverse)
library(ggpubr)

### Load data:
overall_inc_kessler <- vroom("overall_inc_kessler_orig_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
overall_exc_kessler <- vroom("overall_exc_kessler_orig_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
dnmt3a_kessler <- vroom("dnmt3a_kessler_orig_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
tet2_kessler <- vroom("tet2_kessler_orig_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
asxl1_kessler <- vroom("asxl1_kessler_orig_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))

overall_inc_sffdr <- vroom("overall_inc_kessler_sffdr_fp_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
overall_exc_sffdr <- vroom("overall_exc_kessler_sffdr_fp_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
dnmt3a_sffdr <- vroom("dnmt3a_kessler_sffdr_fp_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
tet2_sffdr <- vroom("tet2_kessler_sffdr_fp_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))
asxl1_sffdr <- vroom("asxl1_kessler_sffdr_fp_fuma.txt.gz", col_select = c(SNP, CHR, BP, A1, A2, P))

### Filter SNPs P<0.05 from kessler:
overall_inc_kessler <- overall_inc_kessler %>%
  filter(P<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(overall_inc_kessler)

overall_exc_kessler <- overall_exc_kessler %>%
  filter(P<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(overall_exc_kessler)

dnmt3a_kessler <- dnmt3a_kessler %>%
  filter(P<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(dnmt3a_kessler)

tet2_kessler <- tet2_kessler %>%
  filter(P<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(tet2_kessler)

asxl1_kessler <- asxl1_kessler %>%
  filter(P<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(asxl1_kessler)

### Filter SNPs P<0.05 from sffdr:
overall_inc_sffdr <- overall_inc_sffdr %>%
  rename(FP = P) %>%
  filter(FP<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(overall_inc_sffdr)

overall_exc_sffdr <- overall_exc_sffdr %>%
  rename(FP = P) %>%
  filter(FP<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(overall_exc_sffdr)

dnmt3a_sffdr <- dnmt3a_sffdr %>%
  rename(FP = P) %>%
  filter(FP<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(dnmt3a_sffdr)

tet2_sffdr <- tet2_sffdr %>%
  rename(FP = P) %>%
  filter(FP<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(tet2_sffdr)

asxl1_sffdr <- asxl1_sffdr %>%
  rename(FP = P) %>%
  filter(FP<0.05) %>%
  mutate(id = paste(SNP, CHR, BP, A1, A2, sep = "_")) %>%
  select(-c(SNP, CHR, BP, A1, A2))
dim(asxl1_sffdr)

### Merge datasets:
overall_inc <- overall_inc_kessler %>%
  inner_join(overall_inc_sffdr, by = "id")
dim(overall_inc)

### Define lead and novel variants:
kessler_lead_overall <- c("rs3219104", "rs62175229", "rs6440667", "rs2305407", "rs144317085",
                          "rs7705526", "rs9264393", "rs2396004", "rs2275652", "rs11762008", 
                          "rs7817503", "rs10890839", "rs16930705", "rs72698720", "rs7406485", 
                          "rs80093687", "rs8088824", "rs7272062", "rs2834708", "rs62235635", 
                          "rs361725")
kessler_lead_dnmt3a <- c("rs12745135", "rs2695243", "rs2488000", "rs149934734", "rs61914392", 
                         "rs76428106", "rs2887399", "rs7406485", "rs188761458", "rs8088824", 
                         "rs12472767", "rs6088962", "rs6440667", "rs6773408", "rs144317085", 
                         "rs7705526", "rs2396004", "rs12524502", "rs9389268", "rs11762008",
                         "rs7817503")
kessler_lead_tet2 <- c("rs1505307", "rs2981026", "rs7705526", "rs611646", "rs11846938", "rs78378222")
kessler_lead_asxl1 <- c("rs9374080", "rs2296311")
overall_inc_rsid <- c("rs34016663", "rs1707334", "rs11674811", "rs113542380", "rs342467",
                      "rs1062158", "rs2237154", "rs9404952", "rs536897754", "rs1054684",
                      "rs35150201", "rs55903902", "rs2417347", "rs60626639", "rs4262954",
                      "rs17758695", "rs1474547", "rs1034332")
overall_exc_rsid <- c("rs9267820", "rs12903325")
dnmt3a_rsid <- c("rs1317082", "rs112540634", "rs7800079", "rs11596235", "rs10506957",
                 "rs9535403", "rs7170049", "rs116403919", "rs11668706")
tet2_rsid <- c("rs4918404")
asxl1_rsid <- c("rs4709820", "rs62550973", "rs62227732")

### Categorise novel and known variants:
overall_inc <- overall_inc %>%
  mutate(rsid_extracted = str_extract(id, "^rs\\d+"),
         category = case_when(
           rsid_extracted %in% overall_inc_rsid ~ "Novel",
           rsid_extracted %in% kessler_lead_overall ~ "Kessler lead",
           FP < 5e-8 & P > 5e-8 ~ "Functionally significant", TRUE ~ "Other")) %>%
  select(-rsid_extracted) 

### Calculate -log10():
overall_inc <- overall_inc %>%
  mutate(
    category = factor(category, levels = c("Other", "Functionally significant", "Kessler lead", "Novel")),
    logP = -log10(P),
    logFP = -log10(FP),
    point_size = case_when(category == "Novel" ~ 1.75,
                           category == "Kessler lead" ~ 1.75,
                           category == "Functionally significant" ~ 0.75,
                           category == "Other" ~ 0.75),
    point_alpha = case_when(category == "Other" ~ 0.2,
                            category == "Functionally significant" ~ 0.1, TRUE ~ 1))

### Create plots:
overall_inc_plot <- ggplot() +
  geom_point(data = overall_inc %>% filter(category == "Other"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = overall_inc %>% filter(category == "Functionally significant"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = overall_inc %>% filter(category == "Kessler lead"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = overall_inc %>% filter(category == "Novel"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_vline(xintercept = 7.3, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1.2) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_color_manual(values = c("Novel" = "red2",
                                "Functionally significant" = "plum1",
                                "Kessler lead" = "royalblue1",
                                "Other" = "grey70"),
                     labels = c("Other" = "SNPs",
                                "Functionally significant" = expression(paste("Novel SNPs from sfFDR (FP<5x10"^-8,")")),
                                "Kessler lead" = expression(paste("Lead SNPs from UKBB GWAS (P<5x10"^-8,")")),
                                "Novel" = expression(paste("Novel lead SNPs from sfFDR (FP<5x10"^-8,")")))) + 
  labs(x = expression(-log[10]("P-value from Kessler et al. UKBB GWAS")),
       y = expression(-log[10]("Functional p-value from sfFDR"))) +
  guides(color = guide_legend(title = "SNP category", override.aes = list(size = 5), nrow = 2)) +
  coord_fixed(ratio = 1, xlim = c(0, 20), ylim = c(0, 20), clip = "off") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.justification = c(0.5, 0),
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))
overall_inc_plot
ggsave(overall_inc_plot, file = "kessler_overall_inc_p_vs_fp_plot.pdf")

### Merge datasets:
overall_exc <- overall_exc_kessler %>%
  inner_join(overall_exc_sffdr, by = "id") %>%
  mutate(rsid_extracted = str_extract(id, "^rs\\d+"),
         category = case_when(rsid_extracted %in% overall_exc_rsid ~ "Novel",
                              rsid_extracted %in% kessler_lead_overall ~ "Kessler lead",
                              FP < 5e-8 & P > 5e-8 ~ "Functionally significant", TRUE ~ "Other")) %>%
  select(-rsid_extracted) %>%
  mutate(category = factor(category, levels = c("Other", "Functionally significant", "Kessler lead", "Novel")),
         logP = -log10(P),
         logFP = -log10(FP),
         point_size = case_when(category == "Novel" ~ 1.75,
                                category == "Kessler lead" ~ 1.75,
                                category == "Functionally significant" ~ 0.75,
                                category == "Other" ~ 0.75),
         point_alpha = case_when(category == "Other" ~ 0.2,
                                 category == "Functionally significant" ~ 0.1, TRUE ~ 1))

dnmt3a <- dnmt3a_kessler %>%
  inner_join(dnmt3a_sffdr, by = "id") %>%
  mutate(rsid_extracted = str_extract(id, "^rs\\d+"),
         category = case_when(rsid_extracted %in% dnmt3a_rsid ~ "Novel",
                              rsid_extracted %in% kessler_lead_dnmt3a ~ "Kessler lead",
                              FP < 5e-8 & P > 5e-8 ~ "Functionally significant", TRUE ~ "Other")) %>%
  select(-rsid_extracted) %>%
  mutate(category = factor(category, levels = c("Other", "Functionally significant", "Kessler lead", "Novel")),
         logP = -log10(P),
         logFP = -log10(FP),
         point_size = case_when(category == "Novel" ~ 1.75,
                                category == "Kessler lead" ~ 1.75,
                                category == "Functionally significant" ~ 0.75,
                                category == "Other" ~ 0.75),
         point_alpha = case_when(category == "Other" ~ 0.2,
                                 category == "Functionally significant" ~ 0.1, TRUE ~ 1))

tet2 <- tet2_kessler %>%
  inner_join(tet2_sffdr, by = "id") %>%
  mutate(rsid_extracted = str_extract(id, "^rs\\d+"),
         category = case_when(
           rsid_extracted %in% tet2_rsid ~ "Novel",
           rsid_extracted %in% kessler_lead_tet2 ~ "Kessler lead",
           FP < 5e-8 & P > 5e-8 ~ "Functionally significant",
           TRUE ~ "Other")) %>%
  select(-rsid_extracted) %>%
  mutate(category = factor(category, levels = c("Other", "Functionally significant", "Kessler lead", "Novel")),
         logP = -log10(P),
         logFP = -log10(FP),
         point_size = case_when(category == "Novel" ~ 1.75,
                                category == "Kessler lead" ~ 1.75,
                                category == "Functionally significant" ~ 0.75,
                                category == "Other" ~ 0.75),
         point_alpha = case_when(category == "Other" ~ 0.2,
                                 category == "Functionally significant" ~ 0.1, TRUE ~ 1))

asxl1 <- asxl1_kessler %>%
  inner_join(asxl1_sffdr, by = "id") %>%
  mutate(rsid_extracted = str_extract(id, "^rs\\d+"),
         category = case_when(rsid_extracted %in% asxl1_rsid ~ "Novel",
                              rsid_extracted %in% kessler_lead_asxl1 ~ "Kessler lead",
                              FP < 5e-8 & P > 5e-8 ~ "Functionally significant", TRUE ~ "Other")) %>%
  select(-rsid_extracted) %>%
  mutate(category = factor(category, levels = c("Other", "Functionally significant", "Kessler lead", "Novel")),
         logP = -log10(P),
         logFP = -log10(FP),
         point_size = case_when(category == "Novel" ~ 1.75,
                                category == "Kessler lead" ~ 1.75,
                                category == "Functionally significant" ~ 0.75,
                                category == "Other" ~ 0.75),
         point_alpha = case_when(category == "Other" ~ 0.2,
                                 category == "Functionally significant" ~ 0.1, TRUE ~ 1))

### Create plots and save:
overall_exc_plot <- ggplot() +
  geom_point(data = overall_exc %>% filter(category == "Other"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = overall_exc %>% filter(category == "Functionally significant"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = overall_exc %>% filter(category == "Kessler lead"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = overall_exc %>% filter(category == "Novel"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_vline(xintercept = 7.3, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1.2) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_color_manual(values = c("Novel" = "red2",
                                "Functionally significant" = "plum1",
                                "Kessler lead" = "royalblue1",
                                "Other" = "grey70"),
                     labels = c("Other" = "SNPs",
                                "Functionally significant" = expression(paste("Novel SNPs from sfFDR (FP<5x10"^-8,")")),
                                "Kessler lead" = expression(paste("Lead SNPs from UKBB GWAS (P<5x10"^-8,")")),
                                "Novel" = expression(paste("Novel lead SNPs from sfFDR (FP<5x10"^-8,")")))) + 
  labs(x = expression(-log[10]("P-value from Kessler et al. UKBB GWAS")),
       y = expression(-log[10]("Functional p-value from sfFDR"))) +
  guides(color = guide_legend(title = "SNP category", override.aes = list(size = 5), nrow = 2)) +
  coord_fixed(ratio = 1, xlim = c(0, 20), ylim = c(0, 20), clip = "off") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.justification = c(0.5, 0),
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))
overall_exc_plot
ggsave(overall_exc_plot, file = "kessler_overall_exc_p_vs_fp_plot.pdf")

dnmt3a_plot <- ggplot() +
  geom_point(data = dnmt3a %>%
   filter(category == "Other"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = dnmt3a %>%
   filter(category == "Functionally significant"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = dnmt3a %>%
   filter(category == "Kessler lead"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = dnmt3a %>%
   filter(category == "Novel"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_vline(xintercept = 7.3, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1.2) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_color_manual(values = c("Novel" = "red2",
                                "Functionally significant" = "plum1",
                                "Kessler lead" = "royalblue1",
                                "Other" = "grey70"),
                     labels = c("Other" = "SNPs",
                                "Functionally significant" = expression(paste("Novel SNPs from sfFDR (FP<5x10"^-8,")")),
                                "Kessler lead" = expression(paste("Lead SNPs from UKBB GWAS (P<5x10"^-8,")")),
                                "Novel" = expression(paste("Novel lead SNPs from sfFDR (FP<5x10"^-8,")")))) + 
  labs(x = expression(-log[10]("P-value from Kessler et al. UKBB GWAS")),
       y = expression(-log[10]("Functional p-value from sfFDR"))) +
  guides(color = guide_legend(title = "SNP category", override.aes = list(size = 5), nrow = 2)) +
  coord_fixed(ratio = 1, xlim = c(0, 20), ylim = c(0, 20), clip = "off") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.justification = c(0.5, 0),
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))
dnmt3a_plot
ggsave(dnmt3a_plot, file = "kessler_dnmt3a_p_vs_fp_plot.pdf")

tet2_plot <- ggplot() +
  geom_point(data = tet2 %>%
   filter(category == "Other"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = tet2 %>% 
   filter(category == "Functionally significant"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = tet2 %>% 
   filter(category == "Kessler lead"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = tet2 %>% filter(category == "Novel"),
   aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_vline(xintercept = 7.3, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1.2) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_color_manual(values = c("Novel" = "red2",
                                "Functionally significant" = "plum1",
                                "Kessler lead" = "royalblue1",
                                "Other" = "grey70"),
    labels = c("Other" = "SNPs",
               "Functionally significant" = expression(paste("Novel SNPs from sfFDR (FP<5x10"^-8,")")),
               "Kessler lead" = expression(paste("Lead SNPs from UKBB GWAS (P<5x10"^-8,")")),
               "Novel" = expression(paste("Novel lead SNPs from sfFDR (FP<5x10"^-8,")")))) + 
  labs(x = expression(-log[10]("P-value from Kessler et al. UKBB GWAS")),
       y = expression(-log[10]("Functional p-value from sfFDR"))) +
  guides(color = guide_legend(title = "SNP category", override.aes = list(size = 5), nrow = 2)) +
  coord_fixed(ratio = 1, xlim = c(0, 20), ylim = c(0, 20), clip = "off") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.justification = c(0.5, 0),
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))
tet2_plot
ggsave(tet2_plot, file = "kessler_tet2_p_vs_fp_plot.pdf")

asxl1_plot <- ggplot() +
  geom_point(data = asxl1 %>% filter(category == "Other"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = asxl1 %>% filter(category == "Functionally significant"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = asxl1 %>% filter(category == "Kessler lead"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_point(data = asxl1 %>% filter(category == "Novel"),
             aes(x = logP, y = logFP, size = point_size, alpha = point_alpha, color = category)) +
  geom_vline(xintercept = 7.3, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1.2) +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_color_manual(values = c("Novel" = "red2",
                                "Functionally significant" = "plum1",
                                "Kessler lead" = "royalblue1",
                                "Other" = "grey70"),
                     labels = c("Other" = "SNPs",
                                "Functionally significant" = expression(paste("Novel SNPs from sfFDR (FP<5x10"^-8,")")),
                                "Kessler lead" = expression(paste("Lead SNPs from UKBB GWAS (P<5x10"^-8,")")),
                                "Novel" = expression(paste("Novel lead SNPs from sfFDR (FP<5x10"^-8,")")))) + 
  labs(x = expression(-log[10]("P-value from Kessler et al. UKBB GWAS")),
       y = expression(-log[10]("Functional p-value from sfFDR"))) +
  guides(color = guide_legend(title = "SNP category", override.aes = list(size = 5), nrow = 2)) +
  coord_fixed(ratio = 1, xlim = c(0, 20), ylim = c(0, 20), clip = "off") +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.justification = c(0.5, 0),
        plot.margin = margin(2, 2, 2, 2, unit = "pt"))
asxl1_plot
ggsave(asxl1_plot, file = "kessler_asxl1_p_vs_fp_plot.pdf")