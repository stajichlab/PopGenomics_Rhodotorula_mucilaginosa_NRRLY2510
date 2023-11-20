#!/usr/bin/bash -l 
#SBATCH -p short --mem 2gb --out logs/00_index.log
module load samtools
module load bwa-mem2
if [ -f config.txt ]; then
	source config.txt
fi
mkdir -p $GENOMEFOLDER
pushd $GENOMEFOLDER
# THIS IS EXAMPLE CODE FOR HOW TO DOWNLOAD DIRECT FROM FUNGIDB
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/055/205/GCA_003055205.1_ASM305520v1/GCA_003055205.1_ASM305520v1_genomic.fna.gz

PREF=R_mucilaginosa_JGTA-S1
FASTAFILE=${PREF}_Genome.fasta
## THIS IS FUNGIDB DOWNLOAD PART
echo "working off $FASTAFILE - check if these don't match may need to update config/init script"

if [ ! -f $FASTAFILE ] ; then
	curl -o $FASTAFILE.gz $URL
	pigz -d $FASTAFILE.gz
fi
#if [ ! -f $GFF ]; then
#	curl -O $URL/gff/data/$GFF
#fi

if [[ ! -f $FASTAFILE.fai || $FASTAFILE -nt $FASTAFILE.fai ]]; then
	samtools faidx $FASTAFILE
fi
if [[ ! -f $FASTAFILE.bwt || $FASTAFILE -nt $FASTAFILE.bwt ]]; then
	bwa-mem2 index $FASTAFILE
fi

DICT=$(basename $FASTAFILE .fasta)".dict"

if [[ ! -f $DICT || $FASTAFILE -nt $DICT ]]; then
	rm -f $DICT
	samtools dict $FASTAFILE > $DICT
	ln -s $DICT $FASTAFILE.dict 
fi
#grep ">" $FASTAFILE | perl -p -e 's/>((Chr)?(\d+|mito)_\S+)\s+.+/$1,$3/' > chrom_nums.csv
popd
