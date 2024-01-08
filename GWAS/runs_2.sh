#!/usr/bin/bash -l
#SBATCH -p short -c 32 --mem 96gb 

conda activate ./vc2gwas_env
GFF=/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.genes.csv
echo $GFF
vcf2gwas -v RmucY2510_v1.Phenoytyped.SNP.combined_selected.vcf.gz -pf Rmuc_phenotype_growthrates_regression.vcf2gwas.csv -ap -lmm 4 -gf $GFF -fs 18 -o LMM_run2 -T 31

