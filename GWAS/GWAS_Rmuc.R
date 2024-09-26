#!/usr/bin/env Rscript

library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(tidyverse)
library(dplyr)
library(vcfR)
library(qqman)
library(GAPIT)

bcftoolsexe="/opt/linux/rocky/8.x/x86_64/pkgs/bcftools/1.18/bin/bcftools"
namemap.fn = "../genome/Chrom_Mapping.tab"
vcf.fn <- "RmucY2510_v1.Phenoytyped.SNP.combined_selected.vcf.gz"
vcf.new <- "RmucY2510_v1.Phenoytyped.SNP.good_variants.vcf.gz"
#vcf <- read.vcfR( vcf.fn, verbose = FALSE )
dna_file = "../genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.scaffolds.fa"
dna <- ape::read.dna(dna_file, format = "fasta")

chrom <- create.chromR(name="scaffold", vcf=vcf, seq=dna, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = TRUE)
chrom <- masker(chrom, min_QUAL=0, min_DP=10000, max_DP=25000, min_MQ=59, max_MQ=61)
chrom <- proc.chromR(chrom, verbose = FALSE)
plot(chrom)
write.vcf(chrom, file=vcf.new, mask=TRUE)
MappingNames = read.table(namemap.fn,header=T,stringsAsFactors=FALSE)

pdf(file="GWAS_analysis_Plots.pdf")


pheno <- read_csv("phenotype_growthrates_regression.csv",col_names=TRUE) %>% 
  mutate(scanID = Strain) %>% mutate(sex='F') %>% select(c(scanID,Salt6,T4C,T15C,T22C,T30C,T37C))
# this file needs to be used to create bcftools truncation
outname.good = "RmucY2510_v1.Phenoytyped.SNP.good_variants.renamed.vcf.gz"
print(sprintf("%s annotate --rename-chrs %s %s -Oz -o %s --threads 16",bcftoolsexe,namemap.fn,vcf.new,outname.good))
#system(sprintf("%s annotate --rename-chrs %s %s -Oz -o %s --threads 16",bcftoolsexe,namemap.fn,vcf.new,outname.good))
outname.phenoall = "RmucY2510_v1.Phenoytyped.SNP.combined_selected.renamed.vcf.gz"
print(sprintf("%s annotate --rename-chrs %s %s -Oz -o %s --threads 16",bcftoolsexe,namemap.fn,vcf.fn,outname.phenoall))
#system(sprintf("%s annotate --rename-chrs %s %s -Oz -o %s --threads 16",bcftoolsexe,namemap.fn,vcf.fn,outname.phenoall))
#print(sprintf("bcftools view -Oz -o Rmuc_v7.All.SNP.renamed_truncatedStrains.vcf.gz -S keep_strains_pheno_geno.tsv Rmuc_v7.All.SNP.renamed.vcf.gz


bedfile<-"RmucY2510_v1.Phenoytyped.bed"
bimfile<-"RmucY2510_v1.Phenoytyped.bim"
famfile<-"RmucY2510_v1.Phenoytyped.fam"
gdsfile<-"RmucY2510_v1.Phenoytyped_plink.gds"
gds=snpgdsBED2GDS(bed.fn=bedfile,fam.fn=famfile,bim.fn=bimfile, out.gdsfn=gdsfile,verbose=FALSE)
fam<-read.table(famfile,as.is=TRUE)
names(fam) <-c("family","scanID","father","mother","sex")
fam$sex<-NA
annot<-ScanAnnotationDataFrame(fam)

gdsfile2="RmucY2510_v1.Phenotyped_all.SNP.gds"
# once we already create this
#if (! file.exists(gdsfile))
#  snpgdsVCF2GDS(outname.phenoall, gdsfile,  method="biallelic.only",ignore.chr.prefix = "chr")

snpgdsSummary(gdsfile)

(gds <- GdsGenotypeReader(gdsfile))
strainsphenoname <- pheno$scanID
strainsTyped <- getScanID(gds) 
keepstrains <- tibble(scanID=strainsTyped) %>% inner_join(tibble(scanID=strainsphenoname),by='scanID')
write_tsv(keepstrains,"keep_strains_pheno_geno.tsv")

phenoKeep = data.frame(pheno %>% right_join(keepstrains,by='scanID')) %>% mutate(sex=NA, father=0,mother=0)
phenoScan = ScanAnnotationDataFrame(phenoKeep[order(phenoKeep$scanID),])

getScanID(phenoScan)
close(gds)
sex = rep(NA,length(strainsTyped))
(gds2 <- GdsGenotypeReader(gdsfile2,
                          autosomeCode=1:24L))


snpID <- getSnpID(gds2)
chromosome = as.integer(getChromosome(gds2))
position <- getPosition(gds2)
alleleA <- getAlleleA(gds2)
alleleB <- getAlleleB(gds2)
rsID <- paste(chromosome,position,sep="_")
qual <- getVariable(gds2, "snp.annot/qual")
filter <- getVariable(gds2, "snp.annot/filter")
#qual = rep(30,c(length(snpID)))
filter
qual
snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position,
                                              rsID, alleleA, alleleB,
                                              qual, filter,
                                              stringsAsFactors=FALSE
                                              ))

 
#XchromCode=47L 
#YchromCode=48L
genoData <- GenotypeData(gds2, scanAnnot=phenoScan, snpAnnot=snpAnnot )
getScanVariableNames(genoData)
getGenotype(genoData)


#setVariable(gds,"sample.annot/sex",sex)
#addSex(genoData,sex)
getSex(genoData)

assoc <- assocRegression(genoData, outcome="T37C", model.type="linear")
qqPlot(pval=assoc$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values linear T 37C")
chrom <- getChromosome(snpAnnot, char=TRUE)
snp.sel <- getSnpID(snpAnnot) %in% assoc$snpID
assoc$bp <- with(snpAnnot,position[match(assoc$snpID,snpID)])
manhattanPlot(assoc$Wald.pval, chromosome=chrom[snp.sel],main="Manhattan plot of Wald Test p-values Linear association for T 37C")
manhattan(assoc, chr="chr", bp="bp", snp="snpID", p="Wald.pval" )

# find high scoring SNPs and compare
snp <- pData(snpAnnot)[snp.sel, c("snpID", "rsID")]
snp$pval <- assoc$Wald.pval
snp <- snp[order(snp$pval),]
write.table(snp,"SNP_GWAS_T37C_linear.tsv",sep="\t",quote =FALSE,row.names=FALSE)


assoc <- assocRegression(genoData, outcome="Salt6", model.type="linear")
qqPlot(pval=assoc$Wald.pval, truncate=TRUE, main="QQ Plot of Wald Test p-values linear Salt6")
chrom <- getChromosome(snpAnnot, char=TRUE)
snp.sel <- getSnpID(snpAnnot) %in% assoc$snpID
assoc$bp <- with(snpAnnot,position[match(assoc$snpID,snpID)])
manhattanPlot(assoc$Wald.pval, chromosome=chrom[snp.sel],main="Manhattan plot of Wald Test p-values Linear association for Salt6")
manhattan(assoc, chr="chr", bp="bp", snp="snpID", p="Wald.pval" )

# find high scoring SNPs and compare
snp <- pData(snpAnnot)[snp.sel, c("snpID", "rsID")]
snp$pval <- assoc$Wald.pval
snp <- snp[order(snp$pval),]
write.table(snp,"SNP_GWAS_T37C_linear.tsv",sep="\t",quote =FALSE,row.names=FALSE)


