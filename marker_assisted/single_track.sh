#!/bin/bash

# Requires meryl v1.3 and IGVtools, samtools if fasta is not indexed

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./single_track.sh asm.fasta read.single.meryl k-size"
  exit -1
fi

asm=$1
asm=`echo $asm | sed 's/.gz$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g'`
mer=$2
k=$3

single=$asm.single.k$k

meryl count k=$k $1 output $asm.k$k.meryl
meryl intersect [ equal-to 1 $asm.k$k.meryl ] $mer output $single.meryl
meryl-lookup -bed -sequence $1 -mers $single.meryl > $single.bed
meryl-lookup -wig-depth -sequence $1 -mers $single.meryl > $single.wig

fa=`echo $1 | sed 's/.gz$//g'`
if ! [[ -s $fa.fai ]]; then
  samtools faidx $fa
fi

# Optional, if .tdf is needed
# igvtools count $single.bed $single.tdf $fa.fai

