#!/usr/bin/env Rscript
library(tidyverse)
library(qqman)
tbl = read_tsv("../genome/Chrom_Mapping.tab",col_name=TRUE) %>% mutate(CHR=Old)

phenolist<-c("pH3","pH9","T4C","T15C","T22C","T30C","T37C","NaCL6","NaCL9","NaCL12","Glycerol","LowN","Xylose")
for (PHENO in phenolist) {
results_as <- read.table(sprintf("assoc_results.%s.qassoc", PHENO),head=TRUE)
results = tibble(results_as) %>% left_join(tbl,by="CHR") %>% mutate(CHR=as.integer(New)) %>% select(-c(Old,New)) %>% na.omit()
results
jpeg(sprintf("assoc_manhattan.%s.jpeg",PHENO))
manhattan(data.frame(results),chr="CHR",bp="BP",p="P",snp="SNP", main = sprintf("Manhattan plot: adj assoc %s",PHENO),
	  cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = T, genomewideline = T)
#jpeg(sprintf("QQ-Plot_assoc.%s.jpeg",PHENO))
#qq(results$P, main = sprintf("Q-Q plot of GWAS p-values : log %s",PHENO))
dev.off()
}
