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
	-c	max number of combinatory calls (default 15)
	-o	output prefix
	-y	number of threads
EOF

exit 0

fi

printf "\n\n++++ running: k_eval.sh ++++\n\n"

#set options

while getopts ":a:v:d:k:p:c:o:t:" opt; do

	case $opt in
		a)
			ASM="$OPTARG"
			export ASM
			echo "Assembly: -a $OPTARG"
			;;
		v)
			VAR="$OPTARG"
			export VAR
        		echo "Variants: -v $OPTARG"
			;;
		d)
			DB="$OPTARG"
			echo "Meryl db: -d $OPTARG"
			;;
		k)
			K="$OPTARG"
			export K
			echo "Kmer length: -k $OPTARG"
			;;
		p)
			PEAK="$OPTARG"
			echo "Haploid peak: -p $OPTARG"
			;;
		c)
			COMB="$OPTARG"
			echo "Combinations: -c $OPTARG"
			;;
		o)
			OUT="$OPTARG"
			echo "Output prefix: -o $OPTARG"
			;;   
		t)
			CPU="$OPTARG"
			echo "Number of threads: -t $OPTARG"
			;;         
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;		
	esac
	printf "\n"
done

RAM=$(mktemp -dt "k_eval.XXXXXXXX" --tmpdir=/run/user/$(id -u))
export RAM

printf "Temporary files written to: "
echo $RAM

rm -f ${OUT}.lookup.table

if [[ -z ${COMB} ]]; then

	COMB=15

fi

if [[ ! -e ${OUT}.results ]]; then

	bcftools view --threads $((${CPU}-1)) $VAR -i 'QUAL>20 && (GT="0/1" || GT="1/1")' -Oz > ${RAM}/filtered.vcf.gz
	bcftools index -f ${RAM}/filtered.vcf.gz
	tabix -fp vcf ${RAM}/filtered.vcf.gz

	bcftools view --threads $((${CPU}-1)) -H -i 'QUAL>20 && (GT="0/1" || GT="1/1")' ${RAM}/filtered.vcf.gz > ${RAM}/vars.vcf
	
	n_vars=$(wc -l ${RAM}/vars.vcf | awk '{print $1}')

	a=1
	u=1
	h=1
	STR=1

	while [ $a -le $n_vars ]
	
	do
	
		o=2
	
		truncate -s 0 ${RAM}/var.vcf
		
		var=($(sed -n "${a}p" < ${RAM}/vars.vcf))
		a=$((a+1))

		printf "%s\t"  "${var[@]}" >> ${RAM}/var.vcf
		printf "\n" >> ${RAM}/var.vcf
				
		END=$(( ${var[1]} + ${K} - 2 + ${#var[3]}))
		
		while [[ ! -z "${var}" && $END -gt $STR && $o -le ${COMB} && $a -le $(( $n_vars - 1 )) ]]
		
		do
			
			CHR=${var[0]}
				
			var=($(sed -n "${a}p" < ${RAM}/vars.vcf))
			
			if [[ ! ${var[0]} == $CHR ]]; then
				
				break
				
			fi
			
			STR=$(( ${var[1]} - ${K} + 1))
			
			if [[ $END -gt $STR && $o -le ${COMB} ]]; then

			printf "%s\t"  "${var[@]}" >> ${RAM}/var.vcf
			printf "\n" >> ${RAM}/var.vcf
			
			END=$(( ${var[1]} + ${K} - 2 + ${#var[3]}))		
			
			o=$((o+1))
			a=$((a+1))
			
			fi
					
		done
				
		readarray -t var_ls < ${RAM}/var.vcf
	
		n=${#var_ls[@]}
	
		for (( e = 1; e < 2**n; e++ )); do

			combination=()
	
			idx=$e
	
			for (( j = 0; j < $n; j++ )); do
	
				if (( idx % 2 )); then combination=("${combination[@]}" "${var_ls[$j]}"); fi
	
				idx=$((idx>>1))
	
			done
		
			FIRST=$(echo ${combination[0]} | awk '{print $2}')
			LAST=$(echo ${combination[-1]} | awk '{print $2}')
		
			STR=$(( ${FIRST} - ${K} + 1))
			END=$(( ${LAST} + ${K} - 2 + ${#var[3]}))
			CHR=$(echo ${combination[0]} | awk '{print $1}')

			bcftools view -h $VAR > ${RAM}/var.vcf

			printf "%s\n"  "${combination[@]}" >> ${RAM}/var.vcf
			
			printf "%s\t%s\t%s\t%s\t%s\t%s\n" ${u} ${h} ${CHR} ${STR} ${END} $(bcftools view -H ${RAM}/var.vcf | awk '{printf "%s,",$1":"$2}' | sed 's/.$//' ) >> ${RAM}/${OUT}.lookup.table 
		
			u=$((u+1))

		done
		
		h=$((h+1))

	done
	
	parallel --record-env
	parallel -a ${RAM}/${OUT}.lookup.table -j ${CPU} -k --env _ --colsep '\t' "sh ../k_eval_parallel.sh {1} {3} {4} {5} {6}" >> ${RAM}/var.fa
	
	filename=$(basename -- "${ASM}")
	filename="${filename%.*}"
	
	printf "\n"
	
	meryl count k=${K} ${ASM} output ${filename}.meryl

	meryl-lookup -dump -memory 200 -sequence ${RAM}/var.fa -mers ${filename}.meryl | awk '{print $1"\t"($NF+$(NF-2))}' > ${RAM}/${OUT}.asm
	meryl-lookup -dump -memory 200 -sequence ${RAM}/var.fa -mers ${DB} | awk -v hap=${PEAK} '{COUNT=$NF+$(NF-2); print $1"\t"(COUNT/hap) }' > ${RAM}/${OUT}.read
	
	cut -f2 ${RAM}/${OUT}.read | paste ${RAM}/${OUT}.asm - | awk '{if (index($1, "ref") != 0) {diff=-1} else {diff=1}; 
	if($3==0) {print $0"\tNA\tNA"}
	else if($3>=$2 && $2!=0 && $2+diff!=0){print $0"\t"(-1*($3/$2-1))"\t"(-1*($3/($2+diff)-1))}
	else if($3>=$2 && $2!=0 && $2+diff==0){print $0"\t"(-1*($3/$2-1))"\tNA"}
	else if($3>=$2 && $2==0 && $2+diff!=0){print $0"\tNA\t"(-1*($3/($2+diff)-1))}
	else if($3<$2 && $2!=0 && $2+diff!=0){print $0"\t"($2/$3-1)"\t"(($2+diff)/$3-1)}
	else if($3<$2 && $2!=0 && $2+diff==0){print $0"\t"($2/$3-1)"\tNA"}
	else {printf "whooops"}

	}' > ${RAM}/${OUT}.combined

	awk '{if(NR%2==1) {split($0,var,"_"); if (header!=var[1]) {header=var[1];split(var[1],id,":");split(id[2],coords,"-");printf "%s\t%s\t%s\t%s\t%s\n",substr(header,2),substr(id[1],2),coords[1],coords[2],kmers; kmers=""}}else{kmers=kmers" "$1}}' ${RAM}/var.fa > ${OUT}.kmers.txt

	if [[ -z "$u" ]]; then

		u=$(tail -1 ${RAM}/${OUT}.combined | sed 's/#.*//')
	
	fi

	for ((f=1; f<u; f++))

		do

			ref=($(grep "$(printf "^${f}#")" ${RAM}/${OUT}.combined | grep "ref" | awk 'function abs(x) {return x<0 ? -x : x}
			{
			 if($4!="NA"){pre_sign+=$4;pre_abs_avg+=abs($4)}else{pre_NA+=1};
			 if($5!="NA"){pst_sign+=$5;pst_abs_avg+=abs($5)}else{pst_NA+=1}}
			 END{
			 split($1,a,"#");split(a[2],b,":");split(b[2],c,"_");printf b[1]":"c[1]"\t"
			 if(pre_NA==0){if(pre_sign>=0){var=1}else{var=-1};
			 printf var*pre_abs_avg/NR"\t0"}else{printf "NA\t"pre_NA};
			 if(pst_NA==0){if(pst_sign>=0){var=1}else{var=-1};
			 printf "\t"var*pst_abs_avg/NR"\t0"}else{printf "\tNA\t"pst_NA}}'))
		
			alt=($(grep "$(printf "^${f}#")" ${RAM}/${OUT}.combined | grep "alt" | awk 'function abs(x) {return x<0 ? -x : x}
			{
			 if($4!="NA"){pre_sign+=$4;pre_abs_avg+=abs($4)}else{pre_NA+=1};
			 if($5!="NA"){pst_sign+=$5;pst_abs_avg+=abs($5)}else{pst_NA+=1}}
			 END{
			 if(pre_NA==0){if(pre_sign>=0){var=1}else{var=-1};
			 printf var*pre_abs_avg/NR"\t0"}else{printf "NA\t"pre_NA};
			 if(pst_NA==0){if(pst_sign>=0){var=1}else{var=-1};
			 printf "\t"var*pst_abs_avg/NR"\t0"}else{printf "\tNA\t"pst_NA}}'))

			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${f}" "${ref[0]}" "${ref[1]}" "${ref[2]}" "${ref[3]}" "${ref[4]}" "${alt[0]}" "${alt[1]}" "${alt[2]}" "${alt[3]}" >> ${RAM}/${OUT}.calls

			#awk 'function abs(x) {return x<0 ? -x : x} {if (abs((1-abs($8)))-abs(1-abs($9))>0 && abs((1-abs($10)))-abs(1-abs($11))>0) printf "%s\t%.5f\t%.5f\n", $0,abs((1-abs($8)))-abs(1-abs($9)),abs((1-abs($10)))-abs(1-abs($11))}' ram.calls > ram.diff

	done

	cut -f1,2,6 ${RAM}/${OUT}.lookup.table | paste ${RAM}/${OUT}.calls - > ${OUT}.results

fi