#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                      the T2T project                        ++++
#++++     Credit: Giulio Formenti gformenti@rockefeller.edu       ++++
#++++             Arang  Rhie     arrhie@gmail.com                ++++


bcftools view -Oz ${RAM}/filtered.vcf.gz -r ${5} > ${RAM}/var${1}.vcf.gz
bcftools index -f ${RAM}/var${1}.vcf.gz
tabix -fp vcf ${RAM}/var${1}.vcf.gz

samtools faidx ${ASM} ${2}:${3}-${4} > ${RAM}/ref${1}.fa

REF=($(grep -v ">" ${RAM}/ref${1}.fa | tr -d '\n' | awk -v k=${K} '{for(i=1;i<=length($0)-k+1;i++) print substr($0, 0+i, k)}' ))

for ((i=0; i<"${#REF[@]}"; i++)) do printf ">${1}#${2}:${3}-${4}_ref_%s\n%s\n" $(($i + 1)) "${REF[$i]}"; done

bcftools consensus -f ${RAM}/ref${1}.fa ${RAM}/var${1}.vcf.gz > ${RAM}/alt${1}.fa 2>/dev/null

ALT=($(grep -v ">" ${RAM}/alt${1}.fa | tr -d '\n' | awk -v k=${K} '{for(i=1;i<=length($0)-k+1;i++) print substr($0, 0+i, k)}'))

for ((i=0; i<"${#ALT[@]}"; i++)) do printf ">${1}#${2}:${3}-${4}_alt_%s\n%s\n" $(($i + 1)) "${ALT[$i]}"; done