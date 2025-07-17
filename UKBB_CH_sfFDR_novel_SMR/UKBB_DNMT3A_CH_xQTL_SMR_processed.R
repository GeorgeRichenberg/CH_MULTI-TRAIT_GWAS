### Set working directory and load required packages:
setwd("")
library(vroom)
library(tidyverse)
library(writexl)

### Load data:
dnmt3a_eqtl <- vroom("dnmt3a_trait_eSMR.merged.tsv")
dnmt3a_pqtl <- vroom("dnmt3a_trait_pSMR.merged.tsv")
dnmt3a_mqtl <- vroom("dnmt3a_trait_mSMR.merged.tsv")

### Combine and save:
dnmt3a_mqtl <- dnmt3a_mqtl %>%
  rename(index = Gene)
dnmt3a_xqtl <- rbind(dnmt3a_eqtl, dnmt3a_pqtl, dnmt3a_mqtl)
write_xlsx(dnmt3a_xqtl, "dnmt3a_eqtl_pqtl_mqtl_smr_all.xlsx")

### eQTL:
### Filter significant SMR associations:
max(dnmt3a_eqtl$p_eQTL)
probes_eqtl <- dnmt3a_eqtl %>% distinct(probeID) %>% 
  nrow()
dnmt3a_eqtl2 <- dnmt3a_eqtl %>%
  filter(p_SMR < 1e-3)
probes_eqtl
dim(dnmt3a_eqtl2)

### Filter non-heterogenous associations:
dnmt3a_eqtl3 <- dnmt3a_eqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(dnmt3a_eqtl3)

### pQTL:
### Filter significant SMR associations:
max(dnmt3a_pqtl$p_eQTL)
probes_pqtl <- dnmt3a_pqtl %>% distinct(probeID) %>% 
  nrow()
dnmt3a_pqtl2 <- dnmt3a_pqtl %>%
  filter(p_SMR < 1e-3)
probes_pqtl
dim(dnmt3a_pqtl2)

### Filter non-heterogenous associations:
dnmt3a_pqtl3 <- dnmt3a_pqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(dnmt3a_pqtl3)

### mQTL:
### Filter significant SMR associations:
max(dnmt3a_mqtl$p_eQTL)
probes_mqtl <- dnmt3a_mqtl %>% distinct(probeID) %>% 
  nrow()
dnmt3a_mqtl2 <- dnmt3a_mqtl %>%
  filter(p_SMR < 1e-3)
probes_mqtl
dim(dnmt3a_mqtl2)

### Filter non-heterogenous associations:
dnmt3a_mqtl3 <- dnmt3a_mqtl2 %>%
  filter(p_HEIDI > 0.05)
dim(dnmt3a_mqtl3)

### Save outputs:
write_csv(dnmt3a_eqtl3, file = "dnmt3a_eSMR_sig_non_het.csv")
write_csv(dnmt3a_pqtl3, file = "dnmt3a_pSMR_sig_non_het.csv")
write_csv(dnmt3a_mqtl3, file = "dnmt3a_mSMR_sig_non_het.csv")

### Combine results and save:
dnmt3a_xqtl <- rbind(dnmt3a_eqtl3, dnmt3a_pqtl3, dnmt3a_mqtl3)
write_xlsx(dnmt3a_xqtl, "dnmt3a_eqtl_pqtl_mqtl_smr_sig_1e-3_non_het.xlsx")