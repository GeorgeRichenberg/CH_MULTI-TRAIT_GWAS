### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(BayesRep)
library(purrr)

### Load data:
sffdr_wen <- vroom("sffdr_novel_leadsnps_wen_mcps_compared_bayesrep.csv", col_select = c(CH_GWAS, chr, pos, rsID, nearestGene, beta_Kessler, se_Kessler, beta_Wen, se_Wen))

### Original Kessler effect estimate and standard error:
sffdr_wen <- sffdr_wen %>%
  drop_na()
kessler_beta <- sffdr_wen$beta_Kessler
kessler_se <- sffdr_wen$se_Kessler

### Replication Wen effect estimate and standard error:
wen_beta <- sffdr_wen$beta_Wen
wen_se <- sffdr_wen$se_Wen

### ### Compute and format sceptical Bayes factor and replication Bayes factor for Kessler vs Wen:
sffdr_wen_summary <- sffdr_wen %>% 
  mutate(scepticalBF = BFs(to = kessler_beta, so = kessler_se, tr = wen_beta, sr = wen_se),
         repBF = BFr(to = kessler_beta, so = kessler_se, tr = wen_beta, sr = wen_se)) %>%
  select(CH_GWAS, chr, pos, rsID, nearestGene, beta_Kessler, se_Kessler, beta_Wen, se_Wen, scepticalBF, repBF)
sffdr_wen_summary <- sffdr_wen_summary %>%
  mutate(CH_GWAS = str_remove_all(CH_GWAS, "sfFDR \\[|\\]"))
sffdr_wen_summary

### Save:
write_csv(sffdr_wen_summary, file = "sffdr_novel_kessler_wen_mcps_only_BayesRep.csv")

### Plot posterior distribution of effect size for Kessler vs Wen:
dir.create("BayesRep_kessler_wen_mcps_only_posterior_plots", showWarnings = FALSE)
walk2(seq_along(kessler_beta), sffdr_wen$rsID, ~{png(paste0("BayesRep_kessler_wen_mcps_only_posterior_plots/", .y,
                                                            "_kessler_wen_mcps_only_posterior_plot.png"))
  repPosterior(to = kessler_beta[.x], 
               so = kessler_se[.x], 
               tr = wen_beta[.x], 
               sr = wen_se[.x])
  dev.off()})