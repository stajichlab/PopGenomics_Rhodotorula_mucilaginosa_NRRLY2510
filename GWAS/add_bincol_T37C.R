library(tidyverse)
int <- read_csv("Rmuc_phenotype_growthrates_regression.vcf2gwas.csv")
updated <- int %>% mutate(T37CBinary = case_when(T37C == 0 ~ 0, T37C > 0 ~ 1))
write_csv("Rmuc_phenotype_growthrates_regression.vcf2gwas_bin.csv",updated)
