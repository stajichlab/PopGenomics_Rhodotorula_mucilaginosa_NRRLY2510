#!/usr/bin/env python

import os
import csv
import re
with open("samples.csv",'r') as infh:
    incsv = csv.reader(infh,delimiter=",")
    expected = {}
    i = 1
    for row in incsv:
        if row[0].startswith('Strain') or row[0].startswith('ID'):
            continue
        expected[row[0]] = [i,row[1]]
        i += 1

expectedtmp = expected.copy()
for file in os.listdir("aln"):
    if file.endswith(".cram"):
        pref=file.replace(".cram","")
        if pref in expected:
            expected.pop(pref)
        else:
            print(f"rm -f aln/{file} aln/{file}.crai gvcf/{pref}.g.vcf.gz gvcf/{pref}.g.vcf.gz.tbi")

for file in os.listdir("genomescope"):
    if file.endswith(".histo"):
        pref = file.replace(".histo","")
        if pref not in expectedtmp:
            print(f'rm -rf genomescope/{file} genomescope/{pref}')

for file in os.listdir("gvcf"):
    if file.endswith(".g.vcf.gz"):
        pref = file.replace(".g.vcf.gz","")
        if pref not in expectedtmp:
            print(f'rm -rf gvcf/{file} gvcf/{file}.tbi')


torun = []
indir='input'
base='/bigdata/stajichlab/shared/projects/Rhodotorula/ExtremeRhodotorula_DraftGenomes/input'
for v in expected:
    for fname in expected[v][1].split(';'):
        fname1 = re.sub(r'\[12\]','1',fname)
        if not os.path.exists(os.path.join(indir,fname1)):
#            print(os.path.join(indir,fname1))
            print(f'ln -s {os.path.join(base,fname)} {indir}')
    torun.append(expected[v][0])
torun = ",".join(map( lambda x: str(x), sorted(torun)))
if torun:
    print(f'sbatch -a {torun} pipeline/01_align.sh')
else:
    print('no updates needed')
