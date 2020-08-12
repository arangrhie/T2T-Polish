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
	This script run a titration experiment to evaluate kstar QV at different coverage.
	
	Required inputs are:
	-a	assembly
	-k	kmer length
	-f	forward reads
	-r 	reverse reads
	-1	seq FIRST
	-2	seq INCREMENT
	-3	seq LAST
	-n	number of replicates per condition
	
EOF

exit 0

fi

printf "\n\n++++ running: k_eval.sh ++++\n\n"

#set options

while getopts ":a:k:f:r:1:2:3:n:" opt; do

	case $opt in
		a)
			ASM="$OPTARG"
			echo "Assembly: -a $OPTARG"
			;;
		k)
			K="$OPTARG"
			echo "K length: -k $OPTARG"
			;;
		f)
			FW="$OPTARG"
			echo "Fw reads: -f $OPTARG"
			;;
		r)
			RV="$OPTARG"
        		echo "Rv reads: -r $OPTARG"
			;;
		1)
			FIRST="$OPTARG"
			echo "seq FIRST: -1 $OPTARG"
			;;
		2)
			INCREMENT="$OPTARG"
			echo "seq INCREMENT: -2 $OPTARG"
			;;
		3)
			LAST="$OPTARG"
			echo "seq LAST: -3 $OPTARG"
			;;
		n)
			REP="$OPTARG"
			echo "Number of replicates: -n $OPTARG"
			;;     
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;		
	esac
	printf "\n"
done

mkdir -p logs

meryl count k=${K} ${ASM} output asm.meryl

for exp in $(seq $FIRST $INCREMENT $LAST); do

	for (( rep = 0; rep < $REP; rep++ )); do
	
		log=logs/exp${exp}_rep${rep}.%A_%a.log
		
		if ! [[ -f exp${exp}_rep${rep}.qv ]]; then
	
			sbatch --partition=vgl --cpus-per-task=32 --error=$log --output=$log kstarQV_titration.sh $FW $RV $exp $rep
			
		fi

	done
	
done