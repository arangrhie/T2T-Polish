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
	This script evaluates variant calls based on the K* metric.
	
	Required inputs are:
	-a	assembly multiplicity
	-r 	read multiplicity
	-s	file with chromosome size (assembly.fa.fai)
	-o	output prefix
EOF

exit 0

fi

printf "\n\n++++ running: k_eval.sh ++++\n\n"

#set options

while getopts ":a:v:d:k:p:o:" opt; do

	case $opt in
		a)
			ASM="$OPTARG"
			echo "Assembly: -a $OPTARG"
			;;
		v)
			VAR="$OPTARG"
        		echo "Variants: -v $OPTARG"
			;;
		d)
			DB="$OPTARG"
			echo "Meryl db: -d $OPTARG"
			;;
		k)
			K="$OPTARG"
			echo "Kmer length: -k $OPTARG"
			;;
		p)
			PEAK="$OPTARG"
			echo "Haploid peak: -p $OPTARG"
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

rm -f var.fa

while read -a var

do

if [[ ${#var[3]} == 1 && ${#var[4]} == 1 ]]; then

STR=$(( ${var[1]} - ${K} + 1 ))
END=$(( ${var[1]} + ${K} - 1))

mapfile -t REF < <(vg construct -r ${ASM} -R ${var[0]}:${STR}-${END} | vg kmers -k ${K} - | awk '{split($2,a,":"); if (a[2]>=0 && !(a[2]=="-0")) print a[1]"\t"a[2]"\t"$0}' | sort -nk1 -nk2)
mapfile -t ALT < <(vg construct -r ${ASM} -R ${var[0]}:${STR}-${END} -v ${VAR} | vg kmers -k ${K} - | awk '{split($2,a,":"); if (a[2]>=0 && !(a[2]=="-0")) print a[1]"\t"a[2]"\t"$0}' | sort -nk1 -nk2)

mapfile -t ALT_kmers < <(comm -13 <(printf '%s\n' "${REF[@]}" | awk '{print $3}' | sort) <(printf '%s\n' "${ALT[@]}" | awk '{print $3}' | sort))
printf -v ALT_kmers_joined '%s|' "${ALT_kmers[@]}"
mapfile -t ALT_ONLY < <(printf '%s\n' "${ALT[@]}" | grep -E "${ALT_kmers_joined%|}")

if [[ "${#ALT_ONLY[@]}" -eq "${K}" ]]; then

	printf "\n\nCall:\n\n${var[*]}\n\n"

	printf '\nReference kmers:\n\n'

	for ((idx=0; idx<${#REF[@]}; ++idx))
	do
		e=(${REF[idx]})
		printf "${REF[idx]}\n"
		printf ">${var[0]}:${STR}-${END}_ref_${idx}\n" >> var.fa
		echo ${e[2]} >> var.fa
	done

	printf '\nVariable kmers:\n\n'


	for ((idx=0; idx<${#ALT[@]}; ++idx))
	do
		e=(${ALT[idx]})
		printf "${ALT[idx]}\n"
	done

	printf '\nAlternate kmers:\n\n'

	for ((idx=0; idx<${#ALT_ONLY[@]}; ++idx))
	do
		e=(${ALT_ONLY[idx]})
		printf "${ALT_ONLY[idx]}\n"
		printf ">${var[0]}:${STR}-${END}_alt_${idx}\n" >> var.fa
		echo ${e[2]} >> var.fa
	done
	
fi

fi

done< <(bcftools view -H $VAR)

meryl-lookup -dump -memory 200 -sequence var.fa -mers ../t2t-chm13.20200602.meryl | awk '{print $1"\t"($3+1)"\t"($NF+$(NF-2))}' > ${OUT}.asm
meryl-lookup -dump -memory 200 -sequence var.fa -mers ../CHM13.10X.k21.meryl | awk -v hap=${PEAK} '{ COUNT=$NF+$(NF-2); print $1"\t"($3+1)"\t"(COUNT/hap) }' > ${OUT}.read

cut -f3 ${OUT}.read | paste ${OUT}.asm -