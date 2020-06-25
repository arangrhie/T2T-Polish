#!/bin/bash

in=$1
out=${in/.vcf.gz/}

if [[ -z $in ]]; then
	echo "Usage: ./filt.sh <in.vcf>"
	echo "Requires chr_only.mrg.bed for missing and non_chr_only.bed for non_missing"
	echo "Dependency: java, samtools (bcftools), bedtools"
	exit 0
fi

#q=0.8:alt1
#bcftools view -r chr3,chr8,chr20,chrX -f "PASS" -q $q -Oz 10x_deepvariant.vcf.gz > 10x_deepvariant.chr.PASS.q$q.vcf.gz 
#bcftools view -r chr3,chr8,chr20,chrX -f "PASS" -q $q -Oz hifi_deepvariant.vcf.gz > hifi_deepvariant.chr.PASS.q$q.vcf.gz 
:<<'END'
bcftools view -r chr3,chr8,chr20,chrX -f "PASS" -e 'FORMAT/VAF<0.5 | FORMAT/GQ<30' -Oz 10x_deepvariant.vcf.gz > 10x_deepvariant.PASS.vaf0.5.gq30.vcf.gz 
bcftools view -r chr3,chr8,chr20,chrX -f "PASS" -e 'FORMAT/VAF<0.5 | FORMAT/GQ<30 | FORMAT/DP>70' -Oz 10x_deepvariant.vcf.gz > 10x_deepvariant.PASS.vaf0.5.gq30.dp70.vcf.gz 

bcftools view -R chr_only.mrg.bed -Oz 10x_deepvariant.PASS.vaf0.5.gq30.dp70.vcf.gz > 10x_deepvariant.PASS.vaf0.5.gq30.dp70.missing.vcf.gz
END

bcftools index $in
bcftools view -r chr3,chr8,chr20,chrX -Oz -f "PASS" $in > $out.PASS.vcf.gz
bcftools index $out.PASS.vcf.gz
bcftools view -R chr_only.mrg.bed -Oz $out.PASS.vcf.gz > $out.PASS.missing.vcf.gz
bcftools view -R non_chr_only.bed -Oz $out.PASS.vcf.gz > $out.PASS.non_missing.vcf.gz
java -jar -Xmx1g src/vcfToAlleleCount.jar $out.PASS.missing.vcf.gz $out.PASS.missing.cnt
java -jar -Xmx1g src/vcfToAlleleCount.jar $out.PASS.non_missing.vcf.gz $out.PASS.non_missing.cnt


#module load bedtools
zcat $out.PASS.missing.vcf.gz | grep -v "#" | awk '{if (length($4)>1) {print $1"\t"$2"\t"($2+length($4))} else { print $1"\t"$2"\t"($2+1)}}' | bedtools merge -i - > $out.PASS.missing.bed
zcat $out.PASS.non_missing.vcf.gz | grep -v "#" | awk '{if (length($4)>1) {print $1"\t"$2"\t"($2+length($4))} else { print $1"\t"$2"\t"($2+1)}}' | bedtools merge -i - > $out.PASS.non_missing.bed

ALL=`wc -l chr_only.mrg.bed`
MISSING=`bedtools intersect -a chr_only.mrg.bed -b $out.PASS.missing.bed -u | wc -l`
echo "Variants overlapping missing: $MISSING / $ALL"

