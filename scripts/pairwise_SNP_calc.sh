#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 8 --mem 8gb --out logs/pairwise_SNPs.log

# this script generates number of pairwise differences between strains

module load bcftools
module load workspace/scratch

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
if [ -f config.txt ]; then
	source config.txt
else
	echo "need a config.txt"
fi

if [ -z $FINALVCF ]; then
	echo "Need to define FINALVCF in config.txt"
	exit
fi

if [[ -z $POPYAML || ! -s $POPYAML ]]; then
	echo "Cannot find \$POPYAML variable - set in config.txt"
	exit
fi

OUTDIR=reports/pairwise_strain_compare
mkdir -p $OUTDIR
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
  for TYPE in SNP INDEL
  do
     OUT=$OUTDIR/$PREFIX.$POPNAME.$TYPE.pairwise_count.tsv
     IN=$FINALVCF/$PREFIX.$POPNAME.$TYPE.combined_selected.vcf.gz
    echo -e "STRAIN\tCOUNT" > $OUT
    for strain in $(bcftools query -l $IN)
    do
	    bcftools view --threads $CPU -s $strain -Ob -o $SCRATCH/$strain.bcf $IN
	    bcftools +fill-tags $SCRATCH/$strain.bcf -Ob -o $SCRATCH/$strain.tags.bcf -- -t all
	    count=$(bcftools view -e "AF=0" $SCRATCH/$strain.tags.bcf | grep -c PASS)
	    echo -e "$strain\t$count"
    done >> $OUT
done
