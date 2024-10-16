#!/bin/bash

module load samtools
module load seqtk
module laod aws
module load mashmap # MashMap3
module load bedtools

if [[ $# -lt 3 ]]; then
  echo "Usage: ./collect_summary.sh species asm_version merqury_path"
  echo "Collect assembly summary statistics per haplotypes"
  echo "Final output file will be written as [species]_[asm_version].summary.txt"
  exit -1
fi

sp=$1 # species, e.g. mGorGor1
ver=$2 # assembly version, e.g. v2.0
merqury=$3 # Path to Merqury output; e.g. ./../merqury/${sp_ver}_Hybrid/
sp_ver=${sp}_${ver}
asm=$sp_ver
cpu=$SLURM_CPUS_PER_TASK

if [[ -z $cpu ]]; then
  cpu=8
fi

# Final output: $sp_ver.summary.txt

echo "=== $asm ==="

cat $merqury/*.${sp_ver}.*a*.qv > $sp_ver.qv

echo "Telomere"
seqtk telo -d 10000 $asm.fa > $sp_ver.telo.bed
java -jar -Xmx1g $tools/T2T-Polish/assembly/bedSummarizeTelo.jar $sp_ver.telo.bed $asm.fa.fai > ${sp_ver}.telo.summary

cat ${sp_ver}.telo.summary | grep mat | grep -v "chrM" | \
  grep -v "random" | sort -k1,1V | \
  sed 's/\t3/\tpq/g' | sed 's/\t2/\tq/g' | sed 's/\t1/\tp/g'\
  > $sp_ver.tel.tmp1
cat ${sp_ver}.telo.summary | grep hap1 | grep -v "chrM" | \
  grep -v "random" | sort -k1,1V | \
  sed 's/\t3/\tpq/g' | sed 's/\t2/\tq/g' | sed 's/\t1/\tp/g'\
  >> $sp_ver.tel.tmp1

cat ${sp_ver}.telo.summary | grep pat  | \
  grep -v "random" | sort -k1,1V | awk '{print $2}' | \
  sed 's/3/pq/g' | sed 's/2/q/g' | sed 's/1/p/g' \
  > $sp_ver.tel.tmp2
cat ${sp_ver}.telo.summary | grep hap2 | \
  grep -v "random" | sort -k1,1V | awk '{print $2}' | \
  sed 's/3/pq/g' | sed 's/2/q/g' | sed 's/1/p/g' \
  >> $sp_ver.tel.tmp2

# paste $sp_ver.tel.tmp1 $sp_ver.tel.tmp2

echo "Gap"
seqtk gap $asm.fa > $sp_ver.gap.bed

rm $sp_ver.gap.tmp*
for chr_hap in $(cat $asm.fa.fai | cut -f1 | grep "mat" | grep -v "chrM" | grep -v "random" | sort -k1,1V ) $(cat $asm.fa.fai | cut -f1 | grep "hap1" | grep -v "chrM" | grep -v "random" | sort -k1,1V )
do
  count=`grep -c $chr_hap $sp_ver.gap.bed`
  echo -e "$count" >> $sp_ver.gap.tmp1
done
for chr_hap in $(cat $asm.fa.fai | cut -f1 | grep "pat" | grep -v "random" | sort -k1,1V ) $(cat $asm.fa.fai | cut -f1 | grep "hap2" | grep -v "random" | sort -k1,1V )
do
  count=`grep -c $chr_hap $sp_ver.gap.bed`
  echo -e "$count" >> $sp_ver.gap.tmp2
done
# paste $sp_ver.gap.tmp1 $sp_ver.gap.tmp2

echo "QV"
cat $sp_ver.qv | grep mat  | grep -v "chrM" |\
  grep -v "random" | sort -k1,1V | awk '{print $4}' \
  >  $sp_ver.mqv.tmp1
cat $sp_ver.qv | grep hap1 | grep -v "chrM" |\
  grep -v "random" | sort -k1,1V | awk '{print $4}' \
  >> $sp_ver.mqv.tmp1

cat $sp_ver.qv | grep pat  | \
  grep -v "random" | sort -k1,1V | awk '{print $4}'  > $sp_ver.mqv.tmp2
cat $sp_ver.qv | grep hap2 | \
  grep -v "random" | sort -k1,1V | awk '{print $4}' >> $sp_ver.mqv.tmp2

# paste $sp_ver.mqv.tmp1 $sp_ver.mqv.tmp2

echo
echo "Mashmap 45S"
mashmap \
  -t $cpu \
  --noSplit \
  -q $tools/T2T-Ref/rDNA/human_45S.fa \
  -r $asm.fa \
  -s 13332 \
  --pi 85 \
  -f none \
  -o $sp_ver.45S.mashmap.out

cat $sp_ver.45S.mashmap.out | \
  awk -v OFS='\t' '{print $6, $8, $9, $1, $(NF-1), $5}' | \
  awk -F ":" '{print $1"\t"$NF}' | \
  awk -v OFS='\t' '{print $1, $2, $3, $4, 100*$6, $7}' | \
  awk '$(NF-1)>90' |
  sort -k1,1V -k2,2n > $sp_ver.45S.mashmap.bed

echo
echo "rDNA"
rm $sp_ver.rDNA.tmp*
for chr_hap in $(cat $asm.fa.fai | cut -f1 | grep "mat" | grep -v "chrM" | grep -v "random" | sort -k1,1V ) $(cat $asm.fa.fai | cut -f1 | grep "hap1" | grep -v "chrM" | grep -v "random" | sort -k1,1V )
do
  count=`grep -c $chr_hap $sp_ver.45S.mashmap.bed`
  echo -e "$count" >> $sp_ver.rDNA.tmp1
done
for chr_hap in $(cat $asm.fa.fai | cut -f1 | grep "pat" | grep -v "random" | sort -k1,1V ) $(cat $asm.fa.fai | cut -f1 | grep "hap2" | grep -v "random" | sort -k1,1V )
do
  count=`grep -c $chr_hap $sp_ver.45S.mashmap.bed`
  echo -e "$count" >> $sp_ver.rDNA.tmp2
done
# paste $sp_ver.rDNA.tmp1 $sp_ver.rDNA.tmp2

echo
echo -e "chr\ttel-hap1\ttel-hap2\tgap-hap1\tgap-hap2\trDNA-hap1\trDNA-hap2\tqv-hap1\tqv-hap2"
paste $sp_ver.tel.tmp1 $sp_ver.tel.tmp2 $sp_ver.gap.tmp1 $sp_ver.gap.tmp2 $sp_ver.rDNA.tmp1 $sp_ver.rDNA.tmp2 $sp_ver.mqv.tmp1 $sp_ver.mqv.tmp2 > $sp_ver.summary.txt
cat $sp_ver.summary.txt
echo

cd ../

