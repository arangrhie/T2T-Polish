#!/bin/bash

if [[ "$#" -lt 8 ]]; then
  echo "Usage: ./snv_candidates.sh"
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
if [[ -z $cpus ]]; then
  cpus=12
fi

if [[ -z $SLURM_MEM_PER_NODE ]]; then
  mem="120"
else
  mem=$((SLURM_MEM_PER_NODE/1024))
fi

sequence=$1
seqmer=`echo $sequence | sed 's/\.gz$//' | sed 's/\.fasta$//' | sed 's/\.fa$//'`.meryl

## Create .k31.meryl for sequence
if [[ ! -d  $seqmer ]]; then
  meryl count k=31 $sequence output $seqmer
fi

ver=$2
next_ver=$(echo "$ver" | awk -F "." -v OFS="." '{$NF += 1 ; print}')
echo "Versioning up from $ver to $next_ver"

## Link sequence
if [[ $sequence == *.gz ]]; then
  seq=${ver}.dip.fa.gz
else
  seq=${ver}.dip.fa
fi
mer=$ver.dip.meryl
ln -sf $sequence $seq
ln -sf $seqmer   $mer

## Link vcf files
ont_dip=ont_to_dip
hyb_dip=hybrid_to_dip
hyb_hap1=hybrid_to_hap1
hyb_hap2=hybrid_to_hap2
hyb_hap=hybrid_to_hap

# reheader

module load bcftools

if ! [[ -s $ont_dip.vcf.gz ]]; then
  echo "reheadering vcf files so we can merge them later"
  echo "sample" > name.txt
  bcftools reheader -s name.txt $3 -o $ont_dip.vcf.gz
  bcftools reheader -s name.txt $4 -o $hyb_dip.vcf.gz
  bcftools reheader -s name.txt $5 -o $hyb_hap1.vcf.gz
  bcftools reheader -s name.txt $6 -o $hyb_hap2.vcf.gz
  
  tabix -p vcf $ont_dip.vcf.gz
  tabix -p vcf $hyb_dip.vcf.gz
  tabix -p vcf $hyb_hap1.vcf.gz
  tabix -p vcf $hyb_hap2.vcf.gz
fi

echo "reheader done"

## Link hybrid meryl db of the reads
ln -sf $7 hybrid.meryl
peak=$8


## Filter for PASS
filter_PASS() {
  IN_VCF=$1
  OUT_VCF=${IN_VCF/.vcf.gz/.PASS.vcf.gz}
  bcftools view -f "PASS" \
  -Oz --threads=$cpus  --no-version $IN_VCF > $OUT_VCF
  bcftools index --threads $cpus -f $OUT_VCF
}

set -e
set -o pipefail
set -x

if [[ ! -s $hyb_dip.PASS.vcf.gz ]]; then
  bcftools concat -a -D \
    -Oz --threads=$cpus --no-version \
    $hyb_hap1.vcf.gz $hyb_hap2.vcf.gz > $hyb_hap.vcf.gz
  bcftools index --threads $cpus $hyb_hap.vcf.gz
  
  filter_PASS $hyb_hap.vcf.gz
  filter_PASS $hyb_dip.vcf.gz
  filter_PASS $ont_dip.vcf.gz
  
  ## Consensus error
  bcftools isec -e 'GT="het"' \
    -Oz --threads=$cpus --no-version \
    -p isec_${hyb_hap}.PASS_${hyb_dip}.PASS \
    $hyb_hap.PASS.vcf.gz $hyb_dip.PASS.vcf.gz
  ln -s isec_${hyb_hap}.PASS_${hyb_dip}.PASS/0002.vcf.gz ${hyb_hap}_and_${hyb_dip}.PASS.no_het.vcf.gz
  ln -s isec_${hyb_hap}.PASS_${hyb_dip}.PASS/0002.vcf.gz.tbi \
  ${hyb_hap}_and_${hyb_dip}.PASS.no_het.vcf.gz.tbi
fi

## Merfin
$tools/merfin/1.1/bin/merfin \
  -strict  -peak $peak \
  -sequence $seq \
  -seqmers $mer \
  -memory $mem \
  -output $hyb_dip.PASS.merfin \
  -readmers hybrid.meryl \
  -vcf $hyb_dip.PASS.vcf.gz
mv $hyb_dip.PASS.merfin.filter.vcf $hyb_dip.PASS.merfin-strict.vcf
bcftools view \
  -Oz --threads=$cpus --no-version \
  $hyb_dip.PASS.merfin-strict.vcf > $hyb_dip.PASS.merfin-strict.vcf.gz
bcftools index --threads $cpus $hyb_dip.PASS.merfin-strict.vcf.gz

## Dip homozygous - Hybrid
bcftools view -i 'GT="hom" & GQ>10 & AF>0.8' \
  -Oz --threads=$cpus --no-version \
  $hyb_dip.PASS.vcf.gz > $hyb_dip.PASS.hom_GQ10_AF0.8.vcf.gz
bcftools index --threads $cpus $hyb_dip.PASS.hom_GQ10_AF0.8.vcf.gz

# Hap heterozygous, fixable with ONT support - Hybrid
bcftools view -i 'GT="het" & GQ>20' \
  -Oz --threads=$cpus --no-version \
  $hyb_hap.PASS.vcf.gz > $hyb_hap.PASS.het_GQ20.vcf.gz
bcftools index --threads $cpus $hyb_hap.PASS.het_GQ20.vcf.gz

# Hap heterozygous, fixable with ONT support - ONT
bcftools view -i 'GT="hom" & GQ>20 & AF>0.8' \
  -Oz --threads=$cpus --no-version \
  $ont_dip.PASS.vcf.gz > $ont_dip.PASS.hom_GQ20_AF80.vcf.gz
bcftools index --threads $cpus $ont_dip.PASS.hom_GQ20_AF80.vcf.gz

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

bcftools concat -a -D \
  -Oz --threads=$cpus --no-version \
  ${hyb_hap}_and_${hyb_dip}.PASS.no_het.vcf.gz \
  $hyb_dip.PASS.merfin-strict.vcf.gz \
  $hyb_dip.PASS.hom_GQ10_AF0.8.vcf.gz \
  ${hyb_hap}.het_and_${ont_dip}.hom_and_$hyb_dip.alt.vcf.gz > snv_candidates.vcf.gz
bcftools index --threads $cpus snv_candidates.vcf.gz

# Change GT to 1/1
zcat snv_candidates.vcf.gz | grep "#" > snv_candidates.hom.vcf
zcat snv_candidates.vcf.gz | grep -v "#" | cut -f1-8 | awk '{print $0"\tGT\t1/1"}' >> snv_candidates.hom.vcf

bcftools view -Oz --threads=$cpus --no-version \
  snv_candidates.hom.vcf > snv_candidates.hom.vcf.gz
bcftools index --threads $cpus snv_candidates.hom.vcf.gz

# Merfin -loose
$tools/merfin/1.1/bin/merfin \
  -loose -peak $peak \
  -sequence $seq \
  -seqmers  $mer \
  -memory $mem \
  -output   snv_candidates.merfin \
  -readmers hybrid.meryl \
  -vcf snv_candidates.hom.vcf.gz

# rename
mv snv_candidates.merfin.filter.vcf snv_candidates.merfin-loose.vcf
bcftools view \
  -Oz --threads=$cpus --no-version \
  snv_candidates.merfin-loose.vcf > snv_candidates.merfin-loose.vcf.gz
bcftools index --threads $cpus snv_candidates.merfin-loose.vcf.gz

bcftools consensus -H1 --chain ${ver}_to_${next_ver}.chain \
  -f $seq snv_candidates.merfin-loose.vcf.gz > ${next_ver}.dip.fa
bgzip --index -@$cpus ${next_ver}.dip.fa
samtools faidx -@${cpus} ${next_ver}.dip.fa.gz
