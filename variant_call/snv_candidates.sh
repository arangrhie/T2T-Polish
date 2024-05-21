#!/bin/bash

if [[ "$#" -lt 8 ]]; then
  echo "Usage: ./collect_snv_candidates.sh"
  echo "  sequence.fa \\"
  echo "  ver \\"
  echo "  ont_to_dip.vcf.gz \\"
  echo "  hybrid_to_dip.vcf.gz \\"
  echo "  hybrid_to_hap1.vcf.gz \\"
  echo "  hybrid_to_hap2.vcf.gz \\"
  echo "  hybrid.k31.meryl \\"
  echo "  peak"
  exit -1
fi

cpus=$SLURM_CPUS_PER_TASK

sequence=$1
seqmer=`echo $sequence | sed 's/\.gz$//' | sed 's/\.fasta$//' | sed 's/\.fa$//'`.meryl

## Create .k31.meryl for sequence
if [[ ! -d  $seqmer ]]; then
  meryl count k=31 $sequence output $seqmer
fi

ver=$2

## Link sequence
seq=$ver.analysis-dip.fa
mer=$ver.analysis-dip.meryl
ln -sf $sequence $seq
ln -sf $seqmer   $mer

## Link vcf files
ont_dip=ont_to_dip
hyb_dip=hybrid_to_dip
hyb_hap1=hybrid_to_hap1
hyb_hap2=hybrid_to_hap2
hyb_hap=hybrid_to_hap

ln -sf $3 $ont_dip.vcf.gz
ln -sf $4 $hyb_dip.vcf.gz
ln -sf $5 $hyb_hap1.vcf.gz
ln -sf $6 $hyb_hap2.vcf.gz

ln -sf $3.tbi $ont_dip.vcf.gz.tbi
ln -sf $4.tbi $hyb_dip.vcf.gz.tbi
ln -sf $5.tbi $hyb_hap1.vcf.gz.tbi
ln -sf $6.tbi $hyb_hap2.vcf.gz.tbi

## Link hybrid meryl db of the reads
ln -sf $7 hybrid.meryl
peak=$8

module load bcftools

## Filter for PASS
filter_PASS() {
  IN_VCF=$1
  OUT_VCF=${IN_VCF/.vcf.gz/.PASS.vcf.gz}
  bcftools view -f "PASS" \
  -Oz --threads=$cpus  --no-version $IN_VCF > $OUT_VCF
  bcftools index $OUT_VCF
}

set -e
set -x

bcftools concat -a -D \
  -Oz --threads=$cpus --no-version \
  $hyb_hap1.vcf.gz $hyb_hap2.vcf.gz > $hyb_hap.vcf.gz
bcftools index $hyb_hap.vcf.gz

filter_PASS $hyb_hap.vcf.gz
filter_PASS $hyb_dip.vcf.gz
filter_PASS $ont_dip.vcf.gz

# reheader hybrid_to_dip sample
echo "Hybrid" > hybrid.txt
bcftools reheader -s hybrid.txt $hyb_dip.PASS.vcf.gz > $hyb_dip.PASS.reheader.vcf.gz
mv $hyb_dip.PASS.reheader.vcf.gz $hyb_dip.PASS.vcf.gz
bcftools index $hyb_dip.PASS.vcf.gz

## Consensus error
bcftools isec -e 'GT="het"' \
  -Oz --threads=$cpus --no-version \
  -p isec_${hyb_hap}.PASS_${hyb_dip}.PASS \
  $hyb_hap.PASS.vcf.gz $hyb_dip.PASS.vcf.gz
ln -s isec_${hyb_hap}.PASS_${hyb_dip}.PASS/0002.vcf.gz ${hyb_hap}_and_${hyb_dip}.PASS.no_het.vcf.gz
ln -s isec_${hyb_hap}.PASS_${hyb_dip}.PASS/0002.vcf.gz.tbi \
${hyb_hap}_and_${hyb_dip}.PASS.no_het.vcf.gz.tbi

## Merfin
$tools/merfin/1.1/bin/merfin \
  -strict  -peak $peak \
  -sequence $seq \
  -seqmers $mer \
  -output $hyb_dip.PASS.merfin \
  -readmers hybrid.meryl \
  -vcf $hyb_dip.PASS.vcf.gz
mv $hyb_dip.PASS.merfin.filter.vcf $hyb_dip.PASS.merfin-strict.vcf
bcftools view \
  -Oz --threads=$cpus --no-version \
  $hyb_dip.PASS.merfin-strict.vcf > $hyb_dip.PASS.merfin-strict.vcf.gz
bcftools index $hyb_dip.PASS.merfin-strict.vcf.gz

## Dip homozygous - Hybrid
bcftools view -i 'GT="hom" & GQ>10 & AF>0.8' \
  -Oz --threads=$cpus --no-version \
  $hyb_dip.PASS.vcf.gz > $hyb_dip.PASS.hom_GQ10_AF0.8.vcf.gz
bcftools index $hyb_dip.PASS.hom_GQ10_AF0.8.vcf.gz

# Hap heterozygous, fixable with ONT support - Hybrid
bcftools view -i 'GT="het" & GQ>20' \
  -Oz --threads=$cpus --no-version \
  $hyb_hap.PASS.vcf.gz > $hyb_hap.PASS.het_GQ20.vcf.gz
bcftools index $hyb_hap.PASS.het_GQ20.vcf.gz

# Hap heterozygous, fixable with ONT support - ONT
bcftools view -i 'GT="hom" & GQ>20 & AF>0.8' \
  -Oz --threads=$cpus --no-version \
  $ont_dip.PASS.vcf.gz > $ont_dip.PASS.hom_GQ20_AF80.vcf.gz
bcftools index $ont_dip.PASS.hom_GQ20_AF80.vcf.gz

bcftools isec \
  -p isec_${hyb_hap}.PASS.het_GQ20_${ont_dip}.PASS.hom_GQ20_AF80 \
  -Oz --threads=$cpus --no-version \
  $hyb_hap.PASS.het_GQ20.vcf.gz $ont_dip.PASS.hom_GQ20_AF80.vcf.gz
ln -sf isec_${hyb_hap}.PASS.het_GQ20_${ont_dip}.PASS.hom_GQ20_AF80/0002.vcf.gz ${hyb_hap}.PASS.het_GQ20_and_${ont_dip}.PASS.hom_GQ20_AF80.vcf.gz
ln -sf isec_${hyb_hap}.PASS.het_GQ20_${ont_dip}.PASS.hom_GQ20_AF80/0002.vcf.gz.tbi ${hyb_hap}.PASS.het_GQ20_and_${ont_dip}.PASS.hom_GQ20_AF80.vcf.gz.tbi

# Require one more filter, Hybrid-to-dip "alt"
bcftools view -i 'GT="alt"' \
  -Oz --threads=$cpus --no-version \
  $hyb_dip.PASS.vcf.gz > $hyb_dip.PASS.alt.vcf.gz
bcftools index $hyb_dip.PASS.alt.vcf.gz
bcftools isec \
  -p isec_${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt \
  -Oz --threads=$cpus --no-version \
  ${hyb_hap}.PASS.het_GQ20_and_${ont_dip}.PASS.hom_GQ20_AF80.vcf.gz \
  $hyb_dip.PASS.alt.vcf.gz
ln -sf isec_${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt/0002.vcf.gz ${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt.vcf.gz
ln -sf isec_${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt/0002.vcf.gz.tbi ${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt.vcf.gz.tbi

# Merge to get snv_candidates - hybrid, reheader sample
echo "Sample" > sample.txt
bcftools concat -a -D \
  -Oz --threads=$cpus --no-version \
  ${hyb_hap}_and_${hyb_dip}.PASS.no_het.vcf.gz \
  $hyb_dip.PASS.merfin-strict.vcf.gz \
  $hyb_dip.PASS.hom_GQ10_AF0.8.vcf.gz \
  ${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt.vcf.gz | \
  bcftools reheader -s sample.txt -o snv_candidates.vcf.gz -
bcftools index snv_candidates.vcf.gz

# Change GT to 1/1
zcat snv_candidates.vcf.gz | grep "#" > snv_candidates.hom.vcf
zcat snv_candidates.vcf.gz | grep -v "#" | cut -f1-8 | awk '{print $0"\tGT\t1/1"}' >> snv_candidates.hom.vcf

bcftools view -Oz --threads=$cpus --no-version \
  snv_candidates.hom.vcf > snv_candidates.hom.vcf.gz
bcftools index snv_candidates.hom.vcf.gz

# Merfin -loose
$tools/merfin/1.1/bin/merfin \
  -loose -peak $peak \
  -sequence $seq \
  -seqmers  $mer \
  -output   snv_candidates.merfin \
  -readmers hybrid.meryl \
  -vcf snv_candidates.hom.vcf.gz

# rename
mv snv_candidates.merfin.filter.vcf snv_candidates.merfin-loose.vcf
bcftools view \
  -Oz --threads=$cpus --no-version \
  snv_candidates.merfin-loose.vcf > snv_candidates.merfin-loose.vcf.gz
bcftools index snv_candidates.merfin-loose.vcf.gz
