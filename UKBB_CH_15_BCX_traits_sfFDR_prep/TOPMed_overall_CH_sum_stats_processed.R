### Set working directory and load required pakcages:
setwd("")
library(vroom)
library(tidyverse)
library(data.table)
library(MungeSumstats)

### Load data:
ch1 <- vroom("CHIP_gsr_data_file_chr1_20200604.tsv")
ch2 <- vroom("CHIP_gsr_data_file_chr2_20200604.tsv")
ch3 <- vroom("CHIP_gsr_data_file_chr3_20200604.tsv")
ch4 <- vroom("CHIP_gsr_data_file_chr4_20200604.tsv")
ch5 <- vroom("CHIP_gsr_data_file_chr5_20200604.tsv")
ch6 <- vroom("CHIP_gsr_data_file_chr6_20200604.tsv")
ch7 <- vroom("CHIP_gsr_data_file_chr7_20200604.tsv")
ch8 <- vroom("CHIP_gsr_data_file_chr8_20200604.tsv")
ch9 <- vroom("CHIP_gsr_data_file_chr9_20200604.tsv")
ch10 <- vroom("CHIP_gsr_data_file_chr10_20200604.tsv")
ch11 <- vroom("CHIP_gsr_data_file_chr11_20200604.tsv")
ch12 <- vroom("CHIP_gsr_data_file_chr12_20200604.tsv")
ch13 <- vroom("CHIP_gsr_data_file_chr13_20200604.tsv")
ch14 <- vroom("CHIP_gsr_data_file_chr14_20200604.tsv")
ch15 <- vroom("CHIP_gsr_data_file_chr15_20200604.tsv")
ch16 <- vroom("CHIP_gsr_data_file_chr16_20200604.tsv")
ch17 <- vroom("CHIP_gsr_data_file_chr17_20200604.tsv")
ch18 <- vroom("CHIP_gsr_data_file_chr18_20200604.tsv")
ch19 <- vroom("CHIP_gsr_data_file_chr19_20200604.tsv")
ch20 <- vroom("CHIP_gsr_data_file_chr20_20200604.tsv")
ch21 <- vroom("CHIP_gsr_data_file_chr21_20200604.tsv")
ch22 <- vroom("CHIP_gsr_data_file_chr22_20200604.tsv")
ch23 <- vroom("CHIP_gsr_data_file_chr23_20200604.tsv")

### Combine into one file:
chip_data <- bind_rows(ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9, ch10, 
                       ch11, ch12, ch13, ch14, ch15, ch16, ch17, ch18, ch19, ch20, 
                       ch21, ch22, ch23)
dim(chip_data)
colnames(chip_data)

### Format bick:
chip_data2 <- chip_data %>%
  rename(CHR = Chromosome,
         BP = Position,
         oa_bick = `Allele 1`,
         ea_bick = `Allele 2`, 
         eaf_bick = `Allele 2 frequency`,
         beta_bick = beta,
         se_bick = SE,
         pval_bick = `P value`) %>%
  select(-c(`Coded (effect) allele`, `Sample size`)) %>%
  mutate(CHR = str_remove(CHR, "chr")) %>%
  mutate(CHR = as.numeric(CHR))
chip_data3 <- chip_data2 %>%
  drop_na()
chip_data3 <- filter(chip_data3, eaf_bick >= 0.005)
chip_data3 <- filter(chip_data3, eaf_bick <= 0.995)

### Liftover from GRCh38 to GRCh37 and save:
chip_data3 <- unite(chip_data3, "SNP", CHR:BP, remove = FALSE)
chip_data4 <- data.table(chip_data3)
chip_data_lifted <- liftover(chip_data4, convert_ref_genome = "GRCh37", ref_genome = "GRCh38")
chip_data_lifted <- as_tibble(chip_data_lifted)
dim(chip_data_lifted)
chip_data_lifted <- select(chip_data_lifted, -SNP)
chip_data_lifted <- unite(chip_data_lifted, "SNP", CHR:BP, remove = FALSE)
write_csv(chip_data_lifted, file = "bick_topmed_ch_gwas_buildGRCh37_selectedcols_0.5pcMAF.tsv.gz")
