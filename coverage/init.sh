#!/bin/bash

if [[ -z $1 ]] || [[ -z $2 ]]; then
    echo "Usage: ./init.sh <asm.fasta> <name>"
    echo "  <asm.fasta>: absolute path."
    echo "  <name>: <asm.fasta> will be linked as <name.fa>"
    exit 0
fi

asm=$1
name=$2

module load seqtk

if [[ ! -s $asm.fai ]]; then
  echo "No $asm.fai found. Exit."
  exit -1
fi
echo "# Get sizes"
awk '{print $1"\t0\t"$2}' $asm.fai > $name.bed

echo "# Telomere"
seqtk telo $asm > $name.telo.bed

