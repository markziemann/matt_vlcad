#!/bin/bash
set -x

# Reference transcriptome setup
if [ ! -r ref/gencode.v47.transcripts.fa.idx ] ; then
  mkdir ref
  cd ref
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz
  gunzip gencode.v47.transcripts.fa.gz
  kallisto index -i gencode.v47.transcripts.fa.idx gencode.v47.transcripts.fa
  cd ..
fi

# FASTQC
ls *fastq.gz | parallel fastqc {}

cat Ctrl1*L00*_R1.fastq.gz > Ctl1_R1.fastq.gz
cat Ctrl1*L00*_R2.fastq.gz > Ctl1_R2.fastq.gz

cat Ctrl2*L00*_R1.fastq.gz > Ctl2_R1.fastq.gz
cat Ctrl2*L00*_R2.fastq.gz > Ctl2_R2.fastq.gz

cat Ctrl3*L00*_R1.fastq.gz > Ctl3_R1.fastq.gz
cat Ctrl3*L00*_R2.fastq.gz > Ctl3_R2.fastq.gz

cat Ctrl4*L00*_R1.fastq.gz > Ctl4_R1.fastq.gz
cat Ctrl4*L00*_R2.fastq.gz > Ctl4_R2.fastq.gz

cat Mutant1*L00*_R1.fastq.gz > Mut1_R1.fastq.gz
cat Mutant1*L00*_R2.fastq.gz > Mut1_R2.fastq.gz

cat Mutant2*L00*_R1.fastq.gz > Mut2_R1.fastq.gz
cat Mutant2*L00*_R2.fastq.gz > Mut2_R2.fastq.gz

cat Mutant3*L00*_R1.fastq.gz > Mut3_R1.fastq.gz
cat Mutant3*L00*_R2.fastq.gz > Mut3_R2.fastq.gz

cat Mutant4*L00*_R1.fastq.gz > Mut4_R1.fastq.gz
cat Mutant4*L00*_R2.fastq.gz > Mut4_R2.fastq.gz

IDX=ref/gencode.v47.transcripts.fa.idx

for FQZ1 in Ctl{1..4}*_R1.fastq.gz Mut{1..4}*_R1.fastq.gz ; do
  FQZ2=$(echo $FQZ1 | sed 's#_R1.#_R2.#')
  echo $FQZ1 $FQZ2
  skewer -q 20 -t 16 $FQZ1 $FQZ2
  FQT1=$(echo $FQZ1 | sed 's#fastq.gz#fastq-trimmed-pair1.fastq#')
  FQT2=$(echo $FQZ1 | sed 's#fastq.gz#fastq-trimmed-pair2.fastq#')
  BASE=$(echo $FQZ1 | cut -d '_' -f1)
  kallisto quant -o $BASE -i $IDX -t 16 $FQT1 $FQT2
done

for TSV in $(find . | grep abundance.tsv$) ; do
  NAME=$(echo $TSV | cut -d '/' -f2 )
  cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done | pigz > 3col.tsv.gz

multiqc .
