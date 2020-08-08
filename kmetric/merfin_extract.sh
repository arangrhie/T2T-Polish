#!/bin/bash

bcftools view ${1} -h > ${2%.*}_header.vcf

awk -F "\t" 'function abs(x) {return x<0 ? -x : x}
		{
		 if(id==$2)
		 {
		 if(false<$4){next}
		 else if(false==$4 && abs($7)>abs(trend)){next}		 
		 }
		 else{print var}
		 trend=$7;false=$4;id=$2;var=$9
	}' $2 | sed 's/  /\n/g' | sed '/^$/d' | tr ' ' '\t' > ${2%.*}_vars.vcf

cat ${2%.*}_header.vcf ${2%.*}_vars.vcf > ${2%.*}.merfin.vcf
rm ${2%.*}_header.vcf ${2%.*}_vars.vcf

bcftools view -Oz ${2%.*}.merfin.vcf > ${2%.*}.merfin.vcf.gz
bcftools index ${2%.*}.merfin.vcf.gz
bcftools consensus -f $3 -H 1 ${2%.*}.merfin.vcf.gz > ${3%.*}.polished.fasta
