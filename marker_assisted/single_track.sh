#!/bin/bash

# Requires meryl v1.3 and IGVtools, samtools if fasta is not indexed

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./single_track.sh asm.fasta read.single.meryl k-size"
  exit -1
fi

asm=$1
asm=`echo $asm | sed 's/.fasta$//g' | sed 's/.fa$//g'`
mer=$2
k=$3

meryl count k=$k $1 output $asm.meryl
meryl equal-to 1 $asm.meryl output $asm.1.meryl
meryl intersect $asm.1.meryl $mer output $asm.single.meryl
meryl-lookup -wig-depth -sequence $1 -mers $asm.single.meryl > $asm.single.k$k.wig
meryl-lookup -bed -sequence $1 -mers $asm.single.meryl > $asm.single.k$k.bed

if [[ -z $1.fai ]]; then
  samtools faidx
fi
igvtools count $asm.single.bed $asm.single.tdf $1.fai

