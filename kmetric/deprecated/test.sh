#!/bin/bash

awk 'function abs(x) {return x<0 ? -x : x}
		{
		 if(id==$12)
		 {
		 if(false<$10){next}
		 else if(false==$10 && abs($9)>abs(trend)){next}		 
		 }
		 else{print var}
		 trend=$9;false=$10;id=$12;var=$0
	}' $1 > $1.filtered
	
printf "Predictions\n false kmers lost:%s\n" $(awk '{sum+=$4}END{print sum}' $1.filtered)
	
	
awk '{print $13}' $1.filtered | sed 's/,/\n/g' | sed '1d' | awk '{split($1,a,":");printf "%s\t%s\n",a[1],a[2]}' > $1.vars

bcftools view -h $2 > $1_new.vcf

bcftools view $2 -H | awk 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' $1.vars - >> $1_new.vcf

bcftools view $1_new.vcf -Oz > $1_new.vcf.gz

bcftools index -f $1_new.vcf.gz

bcftools consensus $1_new.vcf.gz -f $3 > ${3#.*}_new.fasta