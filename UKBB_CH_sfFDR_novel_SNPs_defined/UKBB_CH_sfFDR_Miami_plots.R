### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(miamiplot)
library(ggrepel)

### Load data:
overall_inc_sffdr_kessler <- vroom("overall_inc_sffdr_kessler_for_miami.csv", col_types = list(nearestGene = "c"))
overall_exc_sffdr_kessler <- vroom("overall_exc_sffdr_kessler_for_miami.csv", col_types = list(nearestGene = "c"))
dnmt3a_sffdr_kessler <- vroom("dnmt3a_sffdr_kessler_for_miami.csv", col_types = list(nearestGene = "c"))
tet2_sffdr_kessler <- vroom("tet2_sffdr_kessler_for_miami.csv", col_types = list(nearestGene = "c"))
asxl1_sffdr_kessler <- vroom("asxl1_sffdr_kessler_for_miami.csv", col_types = list(nearestGene = "c"))

overall_inc_novel <- vroom("overall_inc_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(rsID, chr, pos, nearestGene))
overall_exc_novel <- vroom("overall_exc_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(rsID, chr, pos, nearestGene))
dnmt3a_novel <- vroom("dnmt3a_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(rsID, chr, pos, nearestGene))
tet2_novel <-vroom("tet2_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(rsID, chr, pos, nearestGene))
asxl1_novel <- vroom("asxl1_sffdr_fuma_novel_leadsnps_vs_kessler_vs_wen_mcps_vs_wen_meta.csv", col_select = c(rsID, chr, pos, nearestGene))

### Re-format as data frames:
overall_inc_sffdr_kessler <- as.data.frame(overall_inc_sffdr_kessler)
overall_exc_sffdr_kessler <- as.data.frame(overall_exc_sffdr_kessler)
dnmt3a_sffdr_kessler <- as.data.frame(dnmt3a_sffdr_kessler)
tet2_sffdr_kessler <- as.data.frame(tet2_sffdr_kessler)
asxl1_sffdr_kessler <- as.data.frame(asxl1_sffdr_kessler)

overall_inc_novel <- as.data.frame(overall_inc_novel)
overall_exc_novel <- as.data.frame(overall_exc_novel)
dnmt3a_novel <- as.data.frame(dnmt3a_novel)
tet2_novel <- as.data.frame(tet2_novel)
asxl1_novel <- as.data.frame(asxl1_novel)

### Rename novel:
overall_inc_novel <- overall_inc_novel %>%
  rename(rsid = rsID)
overall_exc_novel <- overall_exc_novel %>%
  rename(rsid = rsID)
dnmt3a_novel <- dnmt3a_novel %>%
  rename(rsid = rsID)
tet2_novel <- tet2_novel %>%
  rename(rsid = rsID)
asxl1_novel <- asxl1_novel %>%
  rename(rsid = rsID)

### Overall (inc) CH:
### Preprare  data for plotting:
plot_data <- prep_miami_data(data = overall_inc_sffdr_kessler, 
                             split_by = "study", 
                             split_at = "overall_inc_sffdr", 
                             p = "pval")
plot_data_upper <- plot_data$upper %>%
  mutate(is_novel = rsid %in% overall_inc_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - overall_inc_novel$pos[overall_inc_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_upper_axis <- as.data.frame(plot_data_upper %>% 
                                      group_by(chr) %>% 
                                      summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2))
plot_data_upper_novel <- as.data.frame(plot_data_upper %>% 
                                       filter(is_novel == TRUE))
plot_data_upper_near_novel <- as.data.frame(plot_data_upper %>% 
                                            filter(is_near_novel == TRUE))
plot_data_lower <- plot_data$lower %>%
  mutate(is_novel = rsid %in% plot_data_upper_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - overall_inc_novel$pos[overall_inc_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_lower_novel <- plot_data_lower %>%
  filter(rsid %in% plot_data_upper_novel$rsid) %>%
  left_join(plot_data_upper_novel %>% 
  select(rsid, nearestGene), by = "rsid")
plot_data_lower_axis <- as.data.frame(plot_data_lower %>% 
                                      group_by(chr) %>% 
                                      summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2))
plot_data_lower_near_novel <- as.data.frame(plot_data_lower %>% 
                                            filter(is_near_novel == TRUE))

### Create plots:
upper_plot <- ggplot() + 
  geom_point(data = plot_data_upper,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_upper_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) + 
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr, 
                     breaks = plot_data_upper_axis$chr_center,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, max(plot_data_upper$logged_p)),
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('sfFDR analysis: Overall (inc) CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 18),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_upper_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = 5,  
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 10))
upper_plot

lower_plot <- ggplot() + 
  geom_point(data = plot_data_lower,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_lower_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) +
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr,
                     breaks = plot_data_upper_axis$chr_center,
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits  = c(20, 0),
                  expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('UKBB GWAS: Overall (inc) CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 18),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_lower_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene.y),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = -5,
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 10))
lower_plot

### Combine plots and save:
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
p
ggsave(p, file = "overall_inc_sffdr_kessler_plot.pdf")

### Overall (exc) CH:
### Prepare data for plotting:
plot_data <- prep_miami_data(data = overall_exc_sffdr_kessler, 
                             split_by = "study", 
                             split_at = "overall_exc_sffdr", 
                             p = "pval")
plot_data_upper <- plot_data$upper %>%
  mutate(is_novel = rsid %in% overall_exc_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - overall_exc_novel$pos[overall_exc_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_upper_axis <- plot_data_upper %>%
  group_by(chr) %>%
  summarize(chr_center = (max(rel_pos) + min(rel_pos)) / 2) %>%
  as.data.frame()
plot_data_upper_novel <- plot_data_upper %>%
  filter(is_novel == TRUE) %>%
  as.data.frame()
plot_data_upper_near_novel <- plot_data_upper %>%
  filter(is_near_novel == TRUE) %>%
  as.data.frame()
plot_data_lower <- plot_data$lower %>%
  mutate(is_novel = rsid %in% plot_data_upper_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - overall_exc_novel$pos[overall_exc_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_lower_novel <- plot_data_lower %>%
  filter(rsid %in% plot_data_upper_novel$rsid) %>%
  left_join(plot_data_upper_novel %>%
  select(rsid, nearestGene), by = "rsid")
plot_data_lower_near_novel <- plot_data_lower %>%
  filter(is_near_novel == TRUE) %>%
  as.data.frame()

### Create plots:
upper_plot <- ggplot() + 
  geom_point(data = plot_data_upper,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_upper_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) + 
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr, 
                     breaks = plot_data_upper_axis$chr_center,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, max(plot_data_upper$logged_p)), 
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('sfFDR analysis: Overall (exc) CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 18),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_upper_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = 5, 
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 20))

lower_plot <- ggplot() + 
  geom_point(data = plot_data_lower,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_lower_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) +
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr, 
                     breaks = plot_data_upper_axis$chr_center,
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits = c(20, 0),
                  expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('UKBB GWAS: Overall (exc) CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 18),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_lower_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene.y),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = -5,
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 10))
lower_plot

### Combine plots and save:
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
p
ggsave(p, file = "overall_exc_sffdr_kessler_plot.pdf")

### DNMT3A-mutant CH:
### Prepare data for plotting:
plot_data <- prep_miami_data(data = dnmt3a_sffdr_kessler, 
                             split_by = "study", 
                             split_at = "dnmt3a_sffdr", 
                             p = "pval")
plot_data_upper <- plot_data$upper %>%
  mutate(is_novel = rsid %in% dnmt3a_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - dnmt3a_novel$pos[dnmt3a_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_upper_axis <- as.data.frame(plot_data_upper %>% 
                                        group_by(chr) %>% 
                                        summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2))
plot_data_upper_novel <- as.data.frame(plot_data_upper %>% 
                                       filter(is_novel == TRUE))

plot_data_upper_near_novel <- as.data.frame(plot_data_upper %>% 
                                            filter(is_near_novel == TRUE))
plot_data_lower <- plot_data$lower %>%
  mutate(is_novel = rsid %in% plot_data_upper_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - dnmt3a_novel$pos[dnmt3a_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_lower_novel <- plot_data_lower %>%
  filter(rsid %in% plot_data_upper_novel$rsid) %>%
  left_join(plot_data_upper_novel %>%
  select(rsid, nearestGene), by = "rsid")
plot_data_lower_axis <- as.data.frame(plot_data_lower %>% 
                                      group_by(chr) %>% 
                                      summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2))
plot_data_lower_near_novel <- as.data.frame(plot_data_lower %>% 
                                            filter(is_near_novel == TRUE))

### Create plots:
upper_plot <- ggplot() + 
  geom_point(data = plot_data_upper,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_upper_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) + 
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr, 
                     breaks = plot_data_upper_axis$chr_center,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, max(plot_data_upper$logged_p)), 
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('sfFDR analysis: DNMT3A-mutant CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 18),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_upper_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = 5,
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 20))
upper_plot

lower_plot <- ggplot() + 
  geom_point(data = plot_data_lower,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_lower_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) +
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr,
                     breaks = plot_data_upper_axis$chr_center,
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits  = c(20, 0),
                  expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('UKBB GWAS: DNMT3A-mutant CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 18),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_lower_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene.y),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = -5,  
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 10))
lower_plot

### Combine plots and save:
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
p
ggsave(p, file = "dnmt3a_sffdr_kessler_plot.pdf")

### TET2-mutant CH:
### Prepare data for plotting:
plot_data <- prep_miami_data(data = tet2_sffdr_kessler, 
                             split_by = "study", 
                             split_at = "tet2_sffdr", 
                             p = "pval")
plot_data_upper <- plot_data$upper %>%
  mutate(is_novel = rsid %in% tet2_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - tet2_novel$pos[tet2_novel$chr == .x]) <= 500000))) %>%  
  as.data.frame()
plot_data_upper_axis <- plot_data_upper %>% 
  group_by(chr) %>% 
  summarize(chr_center = (max(rel_pos) + min(rel_pos)) / 2) %>% 
  as.data.frame()
plot_data_upper_novel <- plot_data_upper %>% 
  filter(is_novel == TRUE) %>% 
  as.data.frame()
plot_data_upper_near_novel <- plot_data_upper %>% 
  filter(is_near_novel == TRUE) %>% 
  as.data.frame()
plot_data_lower <- plot_data$lower %>%
  mutate(is_novel = rsid %in% plot_data_upper_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - tet2_novel$pos[tet2_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_lower_novel <- plot_data_lower %>%
  filter(rsid %in% plot_data_upper_novel$rsid) %>%
  left_join(plot_data_upper_novel %>%
  select(rsid, nearestGene), by = "rsid")
plot_data_lower_near_novel <- plot_data_lower %>%
  filter(is_near_novel == TRUE) %>%
  as.data.frame()

### Create plots:
upper_plot <- ggplot() + 
  geom_point(data = plot_data_upper,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_upper_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) + 
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr, 
                     breaks = plot_data_lower_axis$chr_center,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, max(plot_data_upper$logged_p)), 
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('sfFDR analysis: TET2-mutant CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 16),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_upper_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = 5,
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 20))
upper_plot

lower_plot <- ggplot() + 
  geom_point(data = plot_data_lower,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_lower_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) +
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr,
                     breaks = plot_data_upper_axis$chr_center,
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits  = c(15, 0),
                 breaks = seq(15, 0, by = -5), 
                  expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('UKBB GWAS: TET2-mutant CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 16),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_lower_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene.y),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = -5,
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold", 
                   ylim = c(NA, 10))
lower_plot

### Combine plots and save:
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
p
ggsave(p, file = "tet2_sffdr_kessler_plot.pdf")

### ASXL1-mutant CH:
### Prepare data for plotting:
plot_data <- prep_miami_data(data = asxl1_sffdr_kessler, 
                             split_by = "study", 
                             split_at = "asxl1_sffdr", 
                             p = "pval")
plot_data_upper <- plot_data$upper %>%
  mutate(is_novel = rsid %in% asxl1_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - asxl1_novel$pos[asxl1_novel$chr == .x]) <= 500000))) %>%  
  as.data.frame()
plot_data_upper_axis <- plot_data_upper %>% 
  group_by(chr) %>% 
  summarize(chr_center = (max(rel_pos) + min(rel_pos)) / 2) %>% 
  as.data.frame()
plot_data_upper_novel <- plot_data_upper %>% 
  filter(is_novel == TRUE) %>% 
  as.data.frame()
plot_data_upper_near_novel <- plot_data_upper %>% 
  filter(is_near_novel == TRUE) %>% 
  as.data.frame()
plot_data_lower <- plot_data$lower %>%
  mutate(is_novel = rsid %in% plot_data_upper_novel$rsid,
         is_near_novel = map2_lgl(chr, pos, ~ any(abs(.y - asxl1_novel$pos[asxl1_novel$chr == .x]) <= 500000))) %>%
  as.data.frame()
plot_data_lower_novel <- plot_data_lower %>%
  filter(rsid %in% plot_data_upper_novel$rsid) %>%
  left_join(plot_data_upper_novel %>% select(rsid, nearestGene), by = "rsid")
plot_data_lower_near_novel <- plot_data_lower %>%
  filter(is_near_novel == TRUE) %>%
  as.data.frame()

### Create plots:
upper_plot <- ggplot() + 
  geom_point(data = plot_data_upper,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_upper_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) + 
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr, 
                     breaks = plot_data_upper_axis$chr_center,
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, ceiling(max(plot_data_upper$logged_p, na.rm = TRUE))), 
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('sfFDR analysis: ASXL1-mutant CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 16),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_upper_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,
                   nudge_y = 1.5,
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold")
upper_plot

lower_plot <- ggplot() + 
  geom_point(data = plot_data_lower,
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)),
             size = 0.75) +
  geom_point(data = plot_data_lower_near_novel,
             aes(x = rel_pos, y = logged_p),
             color = "orange", size = 1) +
  scale_color_manual(values = rep(c("lightblue3", "cornflowerblue"), nrow(plot_data_upper_axis))) + 
  scale_x_continuous(labels = plot_data_upper_axis$chr,
                     breaks = plot_data_upper_axis$chr_center,
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits = c(12, -1),  # Extend lower limit for nudged label
                  breaks = seq(12, 0, by = -3),
                  expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "",
       y = bquote(atop('UKBB GWAS: ASXL1-mutant CH', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),     
        axis.text.x = element_text(size = 14),      
        axis.text.y = element_text(size = 16),  
        plot.margin = margin(b = 0, l = 10)) +
  geom_label_repel(data = plot_data_lower_novel,
                   aes(x = rel_pos, y = logged_p, label = nearestGene.y),
                   size = 5.5, 
                   segment.size = 0.2, 
                   point.padding = 0.3, 
                   force = 2,  
                   nudge_y = -1.5,
                   direction = "y",
                   fill = "orange",
                   box.padding = 0.7,
                   fontface = "bold")
lower_plot

### Combine plots and save:
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
p
ggsave(p, file = "asxl1_sffdr_kessler_plot.pdf")