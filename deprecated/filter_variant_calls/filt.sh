#!/bin/bash

if [[ -z $1 ]]; then
  echo "Usage: filt.sh in.vcf.gz"
  echo "  in.PASS.vcf.gz: PASS filtered variants"
  echo "  in.PASS.cnt: allele counts and other quality values in tdf"
  exit -1
fi

in=$1
out=`echo $in | sed 's/.vcf.gz//' | sed 's/.vcf//'`
out=$out.filt

#q=0.8:alt1
#bcftools view -r chr3,chr8,chr20,chrX -f "PASS" -q $q -Oz 10x_deepvariant.vcf.gz > 10x_deepvariant.chr.PASS.q$q.vcf.gz 

if [[ ! -e $in.tbi ]]; then
  bcftools index $in
fi
bcftools view -q 0.5:alt1 -Oz $in > $out.vcf.gz
java -jar -Xmx1g ~/codes/vcfToAlleleCount.jar $out.vcf.gz > $out.cnt

:<<'END'
#module load bedtools
zcat $out.PASS.missing.vcf.gz | grep -v "#" | awk '{if (length($4)>1) {print $1"\t"$2"\t"($2+length($4))} else { print $1"\t"$2"\t"($2+1)}}' | bedtools merge -i - > $out.PASS.missing.bed
zcat $out.PASS.non_missing.vcf.gz | grep -v "#" | awk '{if (length($4)>1) {print $1"\t"$2"\t"($2+length($4))} else { print $1"\t"$2"\t"($2+1)}}' | bedtools merge -i - > $out.PASS.non_missing.bed

ALL=`wc -l chr_only.mrg.bed`
MISSING=`bedtools intersect -a chr_only.mrg.bed -b $out.PASS.missing.bed -u | wc -l`
NON_MISSING=`bedtools intersect -a chr_only.mrg.bed -b $out.PASS.non_missing.bed -u | wc -l`
echo "Variants overlapping missing: $MISSING / $ALL"
END
