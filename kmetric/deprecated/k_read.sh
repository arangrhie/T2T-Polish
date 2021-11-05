#!/bin/sh

if [[ "$#" -lt 4 ]]; then
	echo "Usage: ./k_read.sh <asm.fasta> <read.meryl> <out> <hap-peak>"
	echo -e "\t<asm.fasta>: assembly fasta file"
	echo -e "\t<read.meryl>: k-mer counts of the reads"
	echo -e "\t<out>: output file prefix. <out>.read will be generated."
	echo -e "\t<hap-peak>: haploid peak multiplicity from read.meryl.histogram"
	exit -1
fi

asm_fa=$1
read=$2
out=$3
hap_peak=$4

echo "
# Dump SEQ-POS SEQ POS COUNT/hap-peak to a table"
# Adjust memory accordingly
meryl-lookup -dump -memory 200 -sequence $asm_fa -mers $read | awk -v hap=$hap_peak '{ COUNT=$NF+$(NF-2); print $1"\t"($3+1)"\t"(COUNT/hap) }' > $out.read

echo "
Done!"
