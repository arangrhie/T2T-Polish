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
	-a	assembly file
	-v 	variant file
	-d	meryl db
	-k	kmer length
	-p	haploid peak
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
rm -f ${OUT}.calls

RAM=$(mktemp -dt "k_eval.XXXXXXXX" --tmpdir=/run/user/$(id -u))

printf "\nTemporary files written to:\n\n"

echo $RAM

while read -a var

do

printf "\n"

printf "%s\t"  "${var[@]}"

printf "\n"

bcftools view -h $VAR > ${RAM}/var.vcf

printf "%s\t"  "${var[@]}" >> ${RAM}/var.vcf

bcftools view -Oz ${RAM}/var.vcf > ${RAM}/var.vcf.gz
bcftools index -f ${RAM}/var.vcf.gz
tabix -fp vcf ${RAM}/var.vcf.gz

STR=$(( ${var[1]} - ${K} + 1 ))
END=$(( ${var[1]} + ${K} - 1))

mapfile -t REF < <(vg construct -r ${ASM} -R ${var[0]}:${STR}-${END} | vg kmers -k ${K} - | awk '{split($2,a,":"); if (a[2]>=0 && !(a[2]=="-0")) print a[1]"\t"a[2]"\t"$0}' | sort -nk1 -nk2)
mapfile -t ALT < <(vg construct -r ${ASM} -R ${var[0]}:${STR}-${END} -v ${RAM}/var.vcf.gz | vg kmers -k ${K} - | awk '{split($2,a,":"); if (a[2]>=0 && !(a[2]=="-0")) print a[1]"\t"a[2]"\t"$0}' | sort -nk1 -nk2)

mapfile -t ALT_kmers < <(comm -13 <(printf '%s\n' "${REF[@]}" | awk '{print $3}' | sort) <(printf '%s\n' "${ALT[@]}" | awk '{print $3}' | sort))
printf -v ALT_kmers_joined '%s|' "${ALT_kmers[@]}"
mapfile -t ALT_ONLY < <(printf '%s\n' "${ALT[@]}" | grep -E "${ALT_kmers_joined%|}")

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
	
done< <(bcftools view -H -i 'QUAL>1 && (GT="0/0" || GT="0/1")' -v snps $VAR)

meryl-lookup -dump -memory 200 -sequence var.fa -mers ../t2t-chm13.20200602.meryl | awk '{print $1"\t"($NF+$(NF-2))}' > ${OUT}.asm
meryl-lookup -dump -memory 200 -sequence var.fa -mers ../CHM13.10X.k21.meryl | awk -v hap=${PEAK} '{ COUNT=$NF+$(NF-2); print $1"\t"(COUNT/hap) }' > ${OUT}.read

cut -f2 ${OUT}.read | paste ${OUT}.asm - | awk '{if (index($1, "ref") != 0) {diff=-1} else {diff=1}; if ($3==0 || $2==0 || $2+1==0 || $2-1==0) {print $0"\t0"} else if ($3>=$2) {print $0"\t"($3/$2)"\t"($3/($2+diff))} else {print $0"\t"(-1*($2/$3))"\t"(-1*(($2+diff)/$3))}}' > ${OUT}.results

while read -a var

do

if [[ ${#var[3]} == 1 && ${#var[4]} == 1 ]]; then

IFS=':' read -r -a arr <<<"${var[9]}"

printf -v call "%s:%s-%s"  "${var[0]}" "$(( ${var[1]} - 20))" "$(( ${var[1]} + 20))"

ref=($(grep "$(echo $call)" ${OUT}.results | grep "ref" | awk '{suma+=$4;sumb+=$5}END{print suma/NR"\t"sumb/NR}'))
alt=($(grep "$(echo $call)" ${OUT}.results | grep "alt" | awk '{suma+=$4;sumb+=$5}END{print suma/NR"\t"sumb/NR}'))

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${var[0]}" "${var[1]}" "${var[3]}"  "${var[4]}"  "${var[5]}" "${arr[0]}" "${arr[1]}" "${ref[0]}" "${ref[1]}" "${alt[0]}" "${alt[1]}" >> ${OUT}.calls

awk 'function abs(x) {return x<0 ? -x : x} {if (abs((1-abs($8)))-abs(1-abs($9))>0 && abs((1-abs($10)))-abs(1-abs($11))>0) printf "%s\t%.5f\t%.5f\n", $0,abs((1-abs($8)))-abs(1-abs($9)),abs((1-abs($10)))-abs(1-abs($11))}' ram.calls > ram.diff

awk '{if(NR%2==1) {split($0,var,"_"); if (header!=var[1]) {header=var[1];split(var[1],id,":");split(id[2],coords,"-");printf "%s\t%s\t%s\t%s\t%s\n",substr(header,2),substr(id[1],2),coords[1],coords[2],kmers; kmers=""}}else{kmers=kmers" "$1}}' var.fa > network.txt
fi

done< <(bcftools view -H -i 'QUAL>1 && (GT="0/0" || GT="0/1")' -v snps $VAR)