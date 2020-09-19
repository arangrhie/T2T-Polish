#!/bin/sh

# Requires samtools, meryl-1.0

asm=$1
samtools faidx $asm.fasta
meryl-lookup -dump -memory 8g -sequence $asm.fasta -mers single.meryl | awk '$4=="T" {print $1"\t"$3"\t"($3+21)}' > $asm.single.bed
igvtools count $asm.single.bed $asm.single.tdf $asm.fasta.fai
