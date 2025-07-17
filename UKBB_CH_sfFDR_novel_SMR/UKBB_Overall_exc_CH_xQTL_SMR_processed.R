### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load data:
overall_exc_eqtl <- vroom("overall_exc_trait_eSMR.merged.tsv")
overall_exc_pqtl <- vroom("overall_exc_trait_pSMR.merged.tsv")
overall_exc_mqtl <- vroom("overall_exc_trait_mSMR.merged.tsv")

### Combine and save:
overall_exc_mqtl <- overall_exc_mqtl %>%
  rename(index = Gene)
overall_exc_xqtl <- rbind(overall_exc_eqtl, overall_exc_pqtl, overall_exc_mqtl)
write_xlsx(overall_exc_xqtl, "overall_exc_eqtl_pqtl_mqtl_smr_all.xlsx")

### eQTL:
### Filter significant SMR associations:
max(overall_exc_eqtl$p_eQTL)
probes_eqtl <- overall_exc_eqtl %>% distinct(probeID) %>% 
  nrow()
overall_exc_eqtl2 <- overall_exc_eqtl %>%
  filter(p_SMR < 1e-3)
probes_eqtl
dim(overall_exc_eqtl2)

### Filter non-heterogenous associations:
overall_exc_eqtl3 <- overall_exc_eqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(overall_exc_eqtl3)

### pQTL:
### Filter significant SMR associations:
max(overall_exc_pqtl$p_eQTL)
probes_pqtl <- overall_exc_pqtl %>% distinct(probeID) %>% 
  nrow()
overall_exc_pqtl2 <- overall_exc_pqtl %>%
  filter(p_SMR < 1e-3)
probes_pqtl
dim(overall_exc_pqtl2)

### Filter non-heterogenous associations:
overall_exc_pqtl3 <- overall_exc_pqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(overall_exc_pqtl3)

### mQTL:
### Filter significant SMR associations:
max(overall_exc_mqtl$p_eQTL)
probes_mqtl <- overall_exc_mqtl %>% distinct(probeID) %>% 
  nrow()
overall_exc_mqtl2 <- overall_exc_mqtl %>%
  filter(p_SMR < 1e-3)
probes_mqtl
dim(overall_exc_mqtl2)

### Filter non-heterogenous associations:
overall_exc_mqtl3 <- overall_exc_mqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(overall_exc_mqtl3)

### Save outputs:
write_csv(overall_exc_eqtl3, file = "overall_exc_eSMR_sig_non_het.csv")
write_csv(overall_exc_pqtl3, file = "overall_exc_pSMR_sig_non_het.csv")
write_csv(overall_exc_mqtl3, file = "overall_exc_mSMR_sig_non_het.csv")

### Combine results and save:
overall_exc_xqtl <- rbind(overall_exc_eqtl3, overall_exc_pqtl3, overall_exc_mqtl3)
write_xlsx(overall_exc_xqtl, "overall_exc_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
