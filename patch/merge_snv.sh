#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./merge_snv.sh species"
  exit -1
fi

module load bedtools
module load samtools

sp=$1

set -x
# new version: v1.4.3
# we need to re-orient and rename the assemblies after this step
OLD=v1.4.1r
NEW=v1.4.3
sp_ver=${sp}_$OLD
sp_new=${sp}_$NEW

# snv edits
SNV=../../variants/$sp/snv_candidates.merfin-loose.vcf.gz
rDNA=../$sp/$sp_ver.rDNA_patch.vcf.gz
TELO=../$sp/$sp_ver.telo_patch.vcf.gz

# Collect regions to exclude
cat ../$sp/target*list | sed 's/,//g' | awk '{print $2}' | awk -F ":" '{print $1"\t"$2"\t"$NF}' | awk -F "-" '{print $1"\t"$NF}' | awk '{print $1"\t"($2-1)"\t"$3}' > $sp_ver.exclude_snv.bed

# Make "Include" bed file
bedtools subtract -a ../../pattern/$sp_ver.bed -b $sp_ver.exclude_snv.bed > $sp_ver.include_snv.bed

bcftools view --no-version -R $sp_ver.include_snv.bed \
  --threads $SLURM_CPUS_PER_TASK -Oz $SNV > $sp_ver.snv_edits.vcf.gz
bcftools index $sp_ver.snv_edits.vcf.gz

if [[ -s $TELO ]]; then
  bcftools concat -aD --no-version \
    --threads $SLURM_CPUS_PER_TASK -Oz $rDNA $TELO $sp_ver.snv_edits.vcf.gz > $sp_ver.snv_sv_edits.vcf.gz
else
  bcftools concat -aD --no-version \
    --threads $SLURM_CPUS_PER_TASK -Oz $rDNA       $sp_ver.snv_edits.vcf.gz > $sp_ver.snv_sv_edits.vcf.gz

fi

bcftools index $sp_ver.snv_sv_edits.vcf.gz

zcat $sp_ver.snv_sv_edits.vcf.gz | grep -v "#" | wc -l

bcftools consensus -c $sp.${OLD}_to_$NEW.chain \
  -f ../$sp/$sp_ver.analysis-dip.fa -HA \
  $sp_ver.snv_sv_edits.vcf.gz > $sp_new.analysis-dip.fa
samtools faidx $sp_new.analysis-dip.fa

# ~/codes/_submit_quick.sh 24 48g $sp_new.qv $MERQURY/eval/qv.sh "../../meryl/${sp}_Hybrid.k31.meryl $sp_new.analysis-dip.fa $sp_new.hybrid"

