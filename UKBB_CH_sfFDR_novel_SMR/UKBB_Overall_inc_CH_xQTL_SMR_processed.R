### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load data:
overall_inc_eqtl <- vroom("overall_inc_trait_eSMR.merged.tsv")
overall_inc_pqtl <- vroom("overall_inc_trait_pSMR.merged.tsv")
overall_inc_mqtl <- vroom("overall_inc_trait_mSMR.merged.tsv")

### Combine and save:
overall_inc_mqtl <- overall_inc_mqtl %>%
  rename(index = Gene)
overall_inc_xqtl <- rbind(overall_inc_eqtl, overall_inc_pqtl, overall_inc_mqtl)
write_xlsx(overall_inc_xqtl, "overall_inc_eqtl_pqtl_mqtl_smr_all.xlsx")

### eQTL:
### Filter significant SMR associations:
max(overall_inc_eqtl$p_eQTL)
probes_eqtl <- overall_inc_eqtl %>% distinct(probeID) %>% 
  nrow()
overall_inc_eqtl2 <- overall_inc_eqtl %>%
  filter(p_SMR < 1e-3)
probes_eqtl
dim(overall_inc_eqtl2)

### Filter non-heterogenous associations:
overall_inc_eqtl3 <- overall_inc_eqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(overall_inc_eqtl3)

### pQTL:
### Filter significant SMR associations:
max(overall_inc_pqtl$p_eQTL)
probes_pqtl <- overall_inc_pqtl %>% distinct(probeID) %>% 
  nrow()
overall_inc_pqtl2 <- overall_inc_pqtl %>%
  filter(p_SMR < 1e-3)
probes_pqtl
dim(overall_inc_pqtl2)

### Filter non-heterogenous associations:
overall_inc_pqtl3 <- overall_inc_pqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(overall_inc_pqtl3)

### mQTL:
### Filter significant SMR associations:
max(overall_inc_mqtl$p_eQTL)
probes_mqtl <- overall_inc_mqtl %>% distinct(probeID) %>% 
  nrow()
overall_inc_mqtl2 <- overall_inc_mqtl %>%
  filter(p_SMR < 1e-3)
probes_mqtl
dim(overall_inc_mqtl2)

### Filter non-heterogenous associations:
overall_inc_mqtl3 <- overall_inc_mqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(overall_inc_mqtl3)

### Save outputs:
write_csv(overall_inc_eqtl3, file = "overall_inc_eSMR_sig_non_het.csv")
write_csv(overall_inc_pqtl3, file = "overall_inc_pSMR_sig_non_het.csv")
write_csv(overall_inc_mqtl3, file = "overall_inc_mSMR_sig_non_het.csv")

### Combine results and save:
overall_inc_xqtl <- rbind(overall_inc_eqtl3, overall_inc_pqtl3, overall_inc_mqtl3)
write_xlsx(overall_inc_xqtl, "overall_inc_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")