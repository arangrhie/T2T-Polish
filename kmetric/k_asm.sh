#!/bin/sh

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./k_asm.sh <asm.fasta> <asm.meryl> <out>"
	echo -e "\t<asm.fasta>: assembly fasta file"
	echo -e "\t<asm.meryl>: assembly meryl dir"
	echo -e "\t<out>: output file prefix. <out>.asm will be generated."
	exit -1
fi

asm_fa=$1
asm=$2
out=$3

echo "
# Dump SEQ-POS SEQ POS COUNT to a table"
meryl-lookup -dump -memory 12 -sequence $asm_fa -mers $asm | awk '$4=="T" {POS=$3+1; print $1"-"POS"\t"$1"\t"POS"\t"($NF+$(NF-2))}' > $out.asm

echo "Done!"

