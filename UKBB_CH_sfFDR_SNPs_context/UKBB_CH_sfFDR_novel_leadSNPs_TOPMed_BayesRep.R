### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(BayesRep)
library(purrr)

### Load data:
sffdr_bick <- vroom("sffdr_novel_leadsnps_bick_compared.csv", col_select = c(CH_GWAS, chr, pos, rsID, nearestGene, beta_kessler, se_kessler, beta_bick, se_bick))

### Original Kessler effect estimate and standard error:
kessler_beta <- sffdr_bick$beta_kessler
kessler_se <- sffdr_bick$se_kessler

### Replication Bick effect estimate and standard error:
bick_beta <- sffdr_bick$beta_bick
bick_se <- sffdr_bick$se_bick

### Compute and format sceptical Bayes factor and replication Bayes factor for Kessler vs Bick:
sffdr_bick_summary <- sffdr_bick %>% 
  mutate(scepticalBF = BFs(to = beta_kessler, so = se_kessler, tr = beta_bick, sr = se_bick),
         repBF = BFr(to = beta_kessler, so = se_kessler, tr = beta_bick, sr = se_bick)) %>%
  select(CH_GWAS, chr, pos, rsID, nearestGene, beta_kessler, se_kessler, beta_bick, se_bick, scepticalBF, repBF)
sffdr_bick_summary <- sffdr_bick_summary %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR \\[|\\]"))
sffdr_bick_summary

### Save:
write_csv(sffdr_bick_summary, file = "sffdr_novel_kessler_bick_BayesRep.csv")

### Plot posterior distribution of effect size for Kessler vs Bick:
dir.create("BayesRep_kessler_bick_posterior_plots", showWarnings = FALSE)
walk2(seq_along(kessler_beta), sffdr_bick$rsID, ~{png(paste0("BayesRep_kessler_bick_posterior_plots/", .y, 
                                                             "_kessler_bick_posterior_plot.png"))
  repPosterior(to = kessler_beta[.x], 
               so = kessler_se[.x], 
               tr = bick_beta[.x], 
               sr = bick_se[.x])
  dev.off()})
