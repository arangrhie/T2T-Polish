#!/bin/bash

if [[ -z $1 ]] || [[ -z $2 ]]; then
    echo "Usage: ./init.sh <asm.fasta> <name>"
    echo "  <asm.fasta>: absolute path."
    echo "  <name>: <asm.fasta> will be linked as <name.fasta>"
    exit 0
fi

asm=$1
name=$2

if [ ! -e $name.fasta ]; then
    echo "Symlink $asm to $name.fasta"
    ln -s $asm $name.fasta
fi

if [ ! -e $name.fasta.fai ]; then
    ln -s $asm.fai $name.fasta.fai
fi

echo "# Get scaffold ends (1kb)"
awk '{print $1"\t0\t"$2}' $name.fasta.fai > $name.bed
cat $name.bed | awk '($3-$2) > 2000 {print $1"\t0\t1000\n"$1"\t"($3-1000)"\t"$3}' > $name.telo.bed  # Ignore scaffolds <2kb

