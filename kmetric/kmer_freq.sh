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
	-c	0-centered (yes|no|both, default yes)
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
		c)
            CENTER="$OPTARG"
			echo "Zero-centered: -c $OPTARG"
			;;          
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;		
	esac
	
printf "\n"

done

if [[ -z  ${CENTER} ]]; then

	CENTER="yes"

fi

#append chr number to generate lookup table

if ! [[ -e "${READ}.table" ]]; then

	awk -F'\t' '{if (!($1 ~ /^[0-9]+$/)) {split($1,chr," ");gsub(/chrom=/, "", chr[2])}else{print chr[2]"-"$1"\t"$0}}' ${READ} > ${READ}.table

fi

if ! [[ -e "${ASM}.table" ]]; then

	awk -F'\t' '{if (!($1 ~ /^[0-9]+$/)) {split($1,chr," ");gsub(/chrom=/, "", chr[2])}else{print chr[2]"-"$1"\t"$0}}' ${ASM} > ${ASM}.table

fi

#join read multiplicity and asm multiplicity in a single table

if ! [[ -e "${OUT}.combined" ]]; then

	join -a 1 ${ASM}.table ${READ}.table | sed 's/ /\t/g' | awk '{if(NF<4) {print $0"\t0\t0"}else{print $0}}' > ${OUT}.combined

fi

if [[  ${CENTER} == "no" ]]; then

	#add autoscaling, calculate exp freq, convert to wig format
	awk -F'\t' -v hap=${PEAK} 'BEGIN{print "track autoScale=on"}{if($2==1){split($1,chrN,"-");print "variableStep chrom="chrN[1]" span=1"};print $2"\t"$5/$3/hap}' ${OUT}.combined > ${OUT}.expfrq.Wig


	#convert to BigWig
	./wigToBigWig ${OUT}.expfrq.Wig ${CHR} ${OUT}.expfrq.bw

elif [ ${CENTER} == "both" ]; then

	#add autoscaling, calculate exp freq, convert to wig format
	awk -F'\t' -v hap=${PEAK} 'BEGIN{print "track autoScale=on"}{if($2==1){split($1,chrN,"-");print "variableStep chrom="chrN[1]" span=1"};print $2"\t"$5/$3/hap}' ${OUT}.combined > ${OUT}.expfrq.Wig

	#generate 0-centered format
	awk '{if (!($1 ~ /^[0-9]+$/)){print $0}else{print $1"\t"$2-1}}' ${OUT}.expfrq.Wig > ${OUT}.expfrq.0.Wig
	
	#convert to BigWig
	wigToBigWig ${OUT}.expfrq.Wig ${CHR} ${OUT}.expfrq.bw
	wigToBigWig ${OUT}.expfrq.0.Wig ${CHR} ${OUT}.expfrq.0.bw

else

	#add autoscaling, calculate exp freq, center on zero, convert to wig format
	awk -F'\t' -v hap=${PEAK} 'BEGIN{print "track autoScale=on"}{if($2==1){split($1,chrN,"-");print "variableStep chrom="chrN[1]" span=1"};print $2"\t"$5/$3/hap-1}' ${OUT}.combined > ${OUT}.expfrq.0.Wig

	#histogram
	awk 'BEGIN{bin_width=0.1}{bin=int(($2-0.0001)/bin_width);if( bin in hist){hist[bin]+=1}else{hist[bin]=1}}END{for (h in hist) printf "%2.2f\t%i\n", h*bin_width, hist[h]}' ${OUT}.expfrq.0.Wig | sort -nk1 > ${OUT}.expfrq.0.Wig.hist

	#convert to BigWig
	wigToBigWig ${OUT}.expfrq.0.Wig ${CHR} ${OUT}.expfrq.0.bw

fi

#upload to globus
#globus transfer XXX:20200602.expfrq.bw YYY:team-curation/kmetric/20200602.expfrq.bw
#globus transfer XXX:20200602.expfrq.bw YYY:team-curation/kmetric/20200602.expfrq.0.bw



