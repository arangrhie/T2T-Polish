#!/bin/sh

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./consensus.sh in vcf out"
	echo -e "\tin\tinput fasta file"
	echo -e "\tvcf\tvcf file containing variants for building the consensus"
	echo -e "\tout\toutput file prefix. out.fasta will be generated"
	exit 0
fi

in=$1
vcf=$2
out=$3

#module load samtools
bcftools consensus -HA -f $in $vcf > $out.fasta
echo "Done. Run Merqury to check remaining missing kmers"

