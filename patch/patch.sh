#!/bin/bash

module load aws/2.15.26
module load seqtk/1.4
module load pangenome/1.1.2

P=500000

# Original assembly to patch
asm_fa=assembly.fa
sp_ver=KOLF2.1Jv0.7  # prefix for the original chr_hap.fa

# Original assembly; since we don't have time and read alignment already exists for this version, let's re-use this
asm_fa=KOLF2.1Jv0.6_maskedY.fa
sp_ver=KOLF2.1Jv0.6

for chr_hap in chr15_2 chr17_1 chr19_2 chr5_2 chr10_1 chr1_2
do
  echo $chr_hap
  # Align $fix_fa to $asm_fa
  if [[ ! -s ${sp_ver}.$chr_hap.fa.fai ]]; then
    samtools faidx $asm_fa $chr_hap > ${sp_ver}.$chr_hap.fa
    samtools faidx ${sp_ver}.$chr_hap.fa
  fi

  out=${chr_hap}_fix_to_${sp_ver}
  echo $out

  wfmash --no-split -ad -t$SLURM_CPUS_PER_TASK ${sp_ver}.$chr_hap.fa \
    ${chr_hap}.fa -s100000 -p95 --one-to-one > $out.sam
  samtools sort -@$SLURM_CPUS_PER_TASK -T $out.tmp -O bam \
    $out.sam > $out.bam
  samtools index $out.bam
  echo
done

# For patching sequences with gaps
## Break patch sequences at gap, collect ~500 kb surrounding the gap
## Link the fix patch assembly as assembly_fix.fa
fix_fa=assembly_fix.fa
fix_gap=assembly_fix.gap.bed

samtools faidx $fix_fa
seqtk gap $fix_fa > $asm_gap
LEN=`wc -l $asm_gap | awk '{print $1}'`

for i in $(seq 1 $LEN)
do
  # Split the fix to contain +- $P sequences surrounding the gap
  chr_hap=`cat $asm_gap | sed -n ${i}p | awk '{print $1}'`
  region=`cat $asm_gap  | sed -n ${i}p | awk -v P=$P '{if ($2-P>0) start=$2-P; else start=1; {print $1":"start"-"$2}}'`
  samtools faidx $fix_fa $region > $chr_hap.1.fa
  region=`cat $asm_gap | sed -n ${i}p | awk -v P=$P '{print $1":"($3+1)"-"($3+1+P)}'`
  samtools faidx $fix_fa $region > $chr_hap.2.fa
  cat $chr_hap.[1-9].fa > ${chr_hap}.both.fa
  cat $chr_hap.[1-2].fa > ${chr_hap}.both.fa
  samtools faidx ${chr_hap}.both.fa
  
  # Align $fix_fa to $asm_fa
  if [[ ! -s samtools faidx ${sp_ver}.$chr_hap.fa.fai ]]; then
    samtools faidx $asm_fa $chr_hap > ${sp_ver}.$chr_hap.fa
    samtools faidx ${sp_ver}.$chr_hap.fa
  fi

  out=${chr_hap}_fix_to_${sp_ver}
  echo $out

  wfmash --no-split -ad -t$SLURM_CPUS_PER_TASK ${sp_ver}.$chr_hap.fa \
    ${chr_hap}.both.fa -s50000 -p95 > $out.sam
  samtools sort -@$SLURM_CPUS_PER_TASK -T $out.tmp -O bam \
    $out.sam > $out.bam
  samtools index $out.bam
  echo
done

samtools merge -@$SLURM_CPUS_PER_TASK -f -o fix_to_${sp_ver}.bam chr*_fix_to_${sp_ver}.bam
samtools index fix_to_${sp_ver}.bam

# Patch alignment ready. Check the boundaries and compare read alignment.
