#!/bin/bash

if [[ "$#" -lt 4 ]]; then
  echo "Usage: run.sh in.bam in.fasta in.single.meryl len_filt [short_reaad]"
  echo "  in.bam: bam file of the whole genome. must be indexed"
  echo "  in.fasta: whole genome fasta"
  echo "  in.single.meryl: marker k-mer db"
  echo "  len_filt: in kb. INTEGER"
  echo "  short_read: Set T if in.bam was generated from short reads. OPTIONAL"
  exit -1
fi

bam=$1
fa=$2
single=$3
len_filt=$4
short_read=$5

if [[ "$short_read" != "T" ]]; then
  echo "Assuming $bam is from long reads"
else
  echo "Assuming $bam is from short reads"
fi

for i in $(seq 1 22) X;
do
echo "
  $tools/T2T-Polish/marker_assisted/_submit.sh $bam chr$i $fa $single $len_filt $short_read;"
  $tools/T2T-Polish/marker_assisted/_submit.sh $bam chr$i $fa $single $len_filt $short_read;
done

