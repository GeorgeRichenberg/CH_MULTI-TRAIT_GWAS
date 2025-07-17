### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load data:
asxl1_eqtl <- vroom("asxl1_trait_eSMR.merged.tsv")
asxl1_pqtl <- vroom("asxl1_trait_pSMR.merged.tsv")
asxl1_mqtl <- vroom("asxl1_trait_mSMR.merged.tsv")

### Combine and save:
asxl1_mqtl <- asxl1_mqtl %>%
  rename(index = Gene)
asxl1_xqtl <- rbind(asxl1_eqtl, asxl1_pqtl, asxl1_mqtl)
write_xlsx(asxl1_xqtl, "asxl1_eqtl_pqtl_mqtl_smr_all.xlsx")

### eQTL:
### Filter significant SMR associations:
max(asxl1_eqtl$p_eQTL)
probes_eqtl <- asxl1_eqtl %>% distinct(probeID) %>% 
  nrow()
asxl1_eqtl2 <- asxl1_eqtl %>%
  filter(p_SMR < 1e-3)
probes_eqtl
dim(asxl1_eqtl2)

### Filter non-heterogenous associations:
asxl1_eqtl3 <- asxl1_eqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(asxl1_eqtl3)

### pQTL:
### Filter significant SMR associations:
max(asxl1_pqtl$p_eQTL)
probes_pqtl <- asxl1_pqtl %>% distinct(probeID) %>% 
  nrow()
asxl1_pqtl2 <- asxl1_pqtl %>%
  filter(p_SMR < 1e-3)
probes_pqtl
dim(asxl1_pqtl2)

### Filter non-heterogenous associations:
asxl1_pqtl3 <- asxl1_pqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(asxl1_pqtl3)

### mQTL:
### Filter significant SMR associations:
max(asxl1_mqtl$p_eQTL)
probes_mqtl <- asxl1_mqtl %>% distinct(probeID) %>% 
  nrow()
asxl1_mqtl2 <- asxl1_mqtl %>%
  filter(p_SMR < 1e-3)
probes_mqtl
dim(asxl1_mqtl2)

### Filter non-heterogenous associations:
asxl1_mqtl3 <- asxl1_mqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(asxl1_mqtl3)

### Save outputs:
write_csv(asxl1_eqtl3, file = "asxl1_eSMR_sig_non_het.csv")
write_csv(asxl1_pqtl3, file = "asxl1_pSMR_sig_non_het.csv")
write_csv(asxl1_mqtl3, file = "asxl1_mSMR_sig_non_het.csv")

### Combine results and save:
asxl1_xqtl <- rbind(asxl1_eqtl3, asxl1_pqtl3, asxl1_mqtl3)
write_xlsx(asxl1_xqtl, "asxl1_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
