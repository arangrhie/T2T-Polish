#!/bin/bash

bcftools annotate -h ${2} ${1} > ${1%.*}.temp.reshaped.vcf

bcftools view -h ${1%.*}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > ${1%.*}.reshaped.vcf

rm ${1%.*}.temp.reshaped.vcf 

bcftools view -H ${1} | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,".\tGT:DP\t1/1:",$8);print $0}' >> ${1%.*}.reshaped.vcf
