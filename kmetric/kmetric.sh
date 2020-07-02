#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                      the T2T project                        ++++
#++++     Credit: Giulio Formenti gformenti@rockefeller.edu       ++++
#++++             Arang  Rhie     arrhie@gmail.com                ++++

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF
	This script generates k* values for all bases in the assembly, given the kmer based copy number estimates of the assembly and in the raw read data.
	
	Required inputs are:
	-a	assembly mutliplicity
	-r 	read multiplicity
	-s	file with chromosome size (assembly.fa.fai)
	-o	output prefix
EOF

exit 0

fi

printf "\n\n++++ running: kmer_freq.sh ++++\n\n"

#set options

while getopts ":a:r:p:s:o:c:" opt; do

	case $opt in
		a)
			ASM="$OPTARG"
			echo "Assembly multiplicity: -a $OPTARG"
			;;
		r)
			READ="$OPTARG"
        		echo "Read multiplicity: -r $OPTARG"
			;;
		s)
			CHR="$OPTARG"
			echo "Chromosome size: -s $OPTARG"
			;;
		o)
			OUT="$OPTARG"
			echo "Output prefix: -o $OPTARG"
			;;           
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;		
	esac
	printf "\n"
done

###join asm and read multiplicity into a single table

# Add header line
echo "track autoScale=on yLineMark=0 yLineOnOff=on" > ${OUT}.expfrq.wig

# READ and ASM are in the same CHR POS COUNT order, so paste them together to CHR POS Ka Kr and get the k* metric
cut -f3 ${READ} | paste ${ASM} - \
	| awk -F '\t' -v CHR="" '{if(CHR!=$1) {print "variableStep chrom="$1" span=1"; CHR=$1}; if ($4==0) {print $2"\t0"} else if ($4>=$3) {print $2"\t"($4/$3)} else {print $2"\t"(-1*($3/$4))}}' \
	>> ${OUT}.expfrq.wig

#convert to BigWig
#module load ucsc
wigToBigWig ${OUT}.expfrq.wig ${CHR} ${OUT}.expfrq.bw

