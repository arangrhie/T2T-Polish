#!/bin/bash

if [[ "$#" -lt 5 ]]; then
  echo "Usage: _submit_marker_filter.sh in.bam in.fasta in.single.meryl len_filt out"
  echo "  in.bam: bam file of the whole genome. must be indexed"
  echo "  in.fasta: whole genome fasta"
  echo "  in.single.meryl: marker k-mer db"
  echo "  len_filt: in kb. INTEGER"
  echo "  out:      output prefix"
  exit -1
fi

bam=$1
fa=$2
single=$3
len_filt=$4
out=$5

for i in $(seq 1 22) X ;
do
  echo "
  sh ~/codes/_submit_norm.sh 24 32g chr$i $tools/T2T-Polish/marker_assisted/filter_by_marker_nosplit.sh  $bam chr$i $fa $single $len_filt"
  sh ~/codes/_submit_norm.sh 24 32g chr$i $tools/T2T-Polish/marker_assisted/filter_by_marker_nosplit.sh "$bam chr$i $fa $single $len_filt" | tail -n1 | awk '{print $NF}' >> filter.jids
done

jids=`cat filter.jids | tr '\n' ',' | sed 's/,$//g'`

echo "
sh ~/codes/_submit_norm.sh 12 24g aggregate $tools/T2T-Polish/marker_assisted/aggregate.sh $out  --dependency=afterok:$jids"
sh ~/codes/_submit_norm.sh 12 24g aggregate $tools/T2T-Polish/marker_assisted/aggregate.sh $out "--dependency=afterok:$jids"
