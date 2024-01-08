#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 24gb 

conda activate ./vcf2gwas_env
#vcf2gwas -v RmucY2510_v1.Phenoytyped.SNP.combined_selected.vcf.gz -pf Rmuc_phenotype_growthrates_regression.vcf2gwas.csv -ap -lmm
vcf2gwas -v RmucY2510_v1.Phenoytyped.SNP.combined_selected.vcf.gz -pf Rmuc_phenotype_growthrates_regression.vcf2gwas.csv -ap -bslmm -fs 18
vcf2gwas -v RmucY2510_v1.Phenoytyped.SNP.combined_selected.vcf.gz -pf Rmuc_phenotype_growthrates_regression.vcf2gwas.csv -ap -lmm -m -fs 18

