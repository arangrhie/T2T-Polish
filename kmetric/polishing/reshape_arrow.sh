#!/bin/bash

bcftools view -h ${1} > ${1%.*}.reshaped.vcf

bcftools view -H ${1} | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,"NULL\tGT:DP\t1/1:",$8);print $0}' >> ${1%.*}.reshaped.vcf
