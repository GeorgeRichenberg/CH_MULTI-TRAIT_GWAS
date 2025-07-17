### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load data:
tet2_eqtl <- vroom("tet2_trait_eSMR.merged.tsv")
tet2_pqtl <- vroom("tet2_trait_pSMR.merged.tsv")
tet2_mqtl <- vroom("tet2_trait_mSMR.merged.tsv")

### Combine and save:
tet2_mqtl <- tet2_mqtl %>%
  rename(index = Gene)
tet2_xqtl <- rbind(tet2_eqtl, tet2_pqtl, tet2_mqtl)
write_xlsx(tet2_xqtl, "tet2_eqtl_pqtl_mqtl_smr_all.xlsx")

### eQTL:
### Filter significant SMR associations:
max(tet2_eqtl$p_eQTL)
probes_eqtl <- tet2_eqtl %>% distinct(probeID) %>% 
  nrow()
tet2_eqtl2 <- tet2_eqtl %>%
  filter(p_SMR < 1e-3)
probes_eqtl
dim(tet2_eqtl2)

### Filter non-heterogenous associations:
tet2_eqtl3 <- tet2_eqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(tet2_eqtl3)

### pQTL:
### Filter significant SMR associations:
max(tet2_pqtl$p_eQTL)
probes_pqtl <- tet2_pqtl %>% distinct(probeID) %>% 
  nrow()
tet2_pqtl2 <- tet2_pqtl %>%
  filter(p_SMR < 1e-3)
probes_pqtl
dim(tet2_pqtl2)

### Filter non-heterogenous associations:
tet2_pqtl3 <- tet2_pqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(tet2_pqtl3)

### mQTL:
### Filter significant SMR associations:
max(tet2_mqtl$p_eQTL)
probes_mqtl <- tet2_mqtl %>% distinct(probeID) %>% 
  nrow()
tet2_mqtl2 <- tet2_mqtl %>%
  filter(p_SMR < 1e-3)
probes_mqtl
dim(tet2_mqtl2)

### Filter non-heterogenous associations:
tet2_mqtl3 <- tet2_mqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(tet2_mqtl3)

### Save outputs:
write_csv(tet2_eqtl3, file = "tet2_eSMR_sig_non_het.csv")
write_csv(tet2_pqtl3, file = "tet2_pSMR_sig_non_het.csv")
write_csv(tet2_mqtl3, file = "tet2_mSMR_sig_non_het.csv")

### Combine results and save:
tet2_xqtl <- rbind(tet2_eqtl3, tet2_pqtl3, tet2_mqtl3)
write_xlsx(tet2_xqtl, "tet2_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")
