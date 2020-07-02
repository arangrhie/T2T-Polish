#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                      the T2T project                        ++++
#++++     Credit: Giulio Formenti gformenti@rockefeller.edu       ++++


if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF
	This script generates kmer equilibrium frequencies for all the bases in the assembly, given the kmer frequency in the assembly and in the raw data.
	
	Inputs are:
	-a	assembly mutliplicity
	-r 	read multiplicity
	-p	haploid peak
	-s	file with chromosome size
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
		p)
            PEAK="$OPTARG"
			echo "Haploid peak: -p $OPTARG"
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

#join read multiplicity and asm multiplicity in a single table

join -a 1 ${ASM} ${READ} | awk '{if(NF<4) {print $0"\t0\t0"}else{print $0}}' | awk -F'\t' -v hap=${PEAK} 'BEGIN{print "track autoScale=on"}{if($3==1){print "variableStep chrom="$2" span=1"};if($5>$3){print $3"\t"$6/$4/hap*-1}else{print $2"\t"$4/$6/hap}' > ${OUT}.expfrq.Wig

#convert to BigWig
./wigToBigWig ${OUT}.expfrq.Wig ${CHR} ${OUT}.expfrq.bw