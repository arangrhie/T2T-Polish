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

rm -f ${OUT}.calls

RAM=$(mktemp -dt "k_eval.XXXXXXXX" --tmpdir=/run/user/$(id -u))

printf "Temporary files written to: "
echo $RAM

if [[ ! -e ${OUT}.results ]]; then

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

		if [[ ! ${var[4]} == *","* ]]; then

			STR=$(( ${var[1]} - ${K} + 1))
			END=$(( ${var[1]} + ${K} - 2 + ${#var[3]}))

			samtools faidx ${ASM} ${var[0]}:${STR}-${END} > ${RAM}/ref.fa

			printf '\nReference sequence:\n\n'

			cat ${RAM}/ref.fa

			printf '\nReference kmers:\n\n'

			REF=($(awk -v k=${K} '{if(NR%2==0) {for(i=1;i<=length($0)-k+1;i++) print substr($0, 0+i, k)}}' ${RAM}/ref.fa))

			paste <(printf "%s\n" "${!REF[@]}" | awk '{print $0+1}') <(printf "%s\n" "${REF[@]}")

			for ((i=0; i< "${#REF[@]}"; i++)) do printf ">${var[0]}:${STR}-${END}_ref_%s\n%s\n" $(($i + 1)) "${REF[$i]}" >> ${RAM}/var.fa; done

			printf '\nAlternate sequence:\n\n'

			sed -E "2s/^(.{$(( ${K} - 1 ))})${var[3]}/\1${var[4]}/" ${RAM}/ref.fa > ${RAM}/alt.fa

			cat ${RAM}/alt.fa

			printf '\nAlternate kmers:\n\n'

			ALT=($(awk -v k=${K} '{if(NR%2==0) {for(i=1;i<=length($0)-k+1;i++) print substr($0, 0+i, k)}}' ${RAM}/alt.fa))

			paste <(printf "%s\n" "${!ALT[@]}" | awk '{print $0+1}') <(printf "%s\n" "${ALT[@]}")

			for ((i=0; i< "${#ALT[@]}"; i++)) do printf ">${var[0]}:${STR}-${END}_alt_%s\n%s\n" $(($i + 1)) "${ALT[$i]}" >> ${RAM}/var.fa; done

		else

			printf '\nComplex variant, skipping.\n\n'

		fi

	done< <(bcftools view -H -i 'QUAL>20 && (GT="0/0" || GT="0/1")' $VAR)

	meryl-lookup -dump -memory 200 -sequence ${RAM}/var.fa -mers ../t2t-chm13.20200602.meryl | awk '{print $1"\t"($NF+$(NF-2))}' > ${OUT}.asm
	meryl-lookup -dump -memory 200 -sequence ${RAM}/var.fa -mers ../CHM13.10X.k21.meryl | awk -v hap=${PEAK} '{ COUNT=$NF+$(NF-2); print $1"\t"(COUNT/hap) }' > ${OUT}.read
	
	cut -f2 ${OUT}.read | paste ${OUT}.asm - | awk '{if (index($1, "ref") != 0) {diff=-1} else {diff=1}; 
	if($3==0) {print $0"\tNA\tNA"}
	else if($3>=$2 && $2!=0 && $2+diff!=0){print $0"\t"(-1*($3/$2-1))"\t"(-1*($3/($2+diff)-1))}
	else if($3>=$2 && $2!=0 && $2+diff==0){print $0"\t"(-1*($3/$2-1))"\tNA"}
	else if($3>=$2 && $2==0 && $2+diff!=0){print $0"\tNA\t"(-1*($3/($2+diff)-1))}
	else if($3<$2 && $2!=0 && $2+diff!=0){print $0"\t"($2/$3-1)"\t"(($2+diff)/$3-1)}
	else if($3<$2 && $2!=0 && $2+diff==0){print $0"\t"($2/$3-1)"\tNA"}
	else if($3<$2 && $2==0 && $2+diff!=0){print $0"\tNA\t"(($2+diff)/$3-1)}

	}' > ${OUT}.results

	awk '{if(NR%2==1) {split($0,var,"_"); if (header!=var[1]) {header=var[1];split(var[1],id,":");split(id[2],coords,"-");printf "%s\t%s\t%s\t%s\t%s\n",substr(header,2),substr(id[1],2),coords[1],coords[2],kmers; kmers=""}}else{kmers=kmers" "$1}}' ${RAM}/var.fa > network.txt

fi

while read -a var

	do

	if [[ ! ${var[4]} == *","* ]]; then

		IFS=':' read -r -a arr <<<"${var[9]}"

		printf -v call "%s:%s-%s"  "${var[0]}" "$(( ${var[1]} - ${K} + 1))" "$(( ${var[1]} + ${K} - 2 + ${#var[3]}))"

		ref=($(grep "$(echo $call)" ${OUT}.results | grep "ref" | awk 'function abs(x) {return x<0 ? -x : x} {if($4!="NA"){suma+=$4;avga+=abs($4)}else{suma="NA"};if($5!="NA"){sumb+=$5;avgb+=abs($5)}else{sumb="NA"}}END{if(suma!="NA"){if(suma>0){var=1}else{var=-1};printf var*avga/NR}else{printf "NA"};if(sumb!="NA"){if(sumb>0){var=1}else{var=-1};printf "\t"var*avgb/NR}else{printf "\tNA"}}'))
		alt=($(grep "$(echo $call)" ${OUT}.results | grep "alt" | awk 'function abs(x) {return x<0 ? -x : x} {if($4!="NA"){suma+=$4;avga+=abs($4)}else{suma="NA"};if($5!="NA"){sumb+=$5;avgb+=abs($5)}else{sumb="NA"}}END{if(suma!="NA"){if(suma>0){var=1}else{var=-1};printf var*avga/NR}else{printf "NA"};if(sumb!="NA"){if(sumb>0){var=1}else{var=-1};printf "\t"var*avgb/NR}else{printf "\tNA"}}'))

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${var[0]}" "${var[1]}" "${var[3]}"  "${var[4]}"  "${var[5]}" "${arr[0]}" "${arr[1]}" "${ref[0]}" "${ref[1]}" "${alt[0]}" "${alt[1]}" >> ${OUT}.calls

		awk 'function abs(x) {return x<0 ? -x : x} {if (abs((1-abs($8)))-abs(1-abs($9))>0 && abs((1-abs($10)))-abs(1-abs($11))>0) printf "%s\t%.5f\t%.5f\n", $0,abs((1-abs($8)))-abs(1-abs($9)),abs((1-abs($10)))-abs(1-abs($11))}' ram.calls > ram.diff

	else

		printf '\nComplex variant, skipping.\n\n'

	fi

done< <(bcftools view -H -i 'QUAL>20 && (GT="0/0" || GT="0/1")' $VAR)