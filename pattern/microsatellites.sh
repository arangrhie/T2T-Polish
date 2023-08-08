#!/bin/bash

if [[ -z $1 ]]; then
  echo "Usage: ./microsatellites.sh in.fasta"
  echo "Collect microsatellite repeat patterns. Output files will be written in <in> dir."
  exit -1
fi
fa=$1
out=`echo $fa | sed 's/\.fasta$//g' | sed 's/\.fa$//g' | sed 's/\.fasta\.gz$//g' | sed 's/\.fa\.gz$//g'`

module load bedtools

mkdir -p $out

for p in ga gc at
do
  echo "Collect $p"
  $tools/seqrequester/build/bin/seqrequester microsatellite -prefix $out/$out.microsatellite -window 128 -$p $fa
done

for p in GC AT GA TC
do
  echo "Convert to bed: $p"
  # java -jar -Xmx1g $tools/T2T-Polish/paf_util/bedToFixedWig.jar $out/$out.microsatellite.$p.128.bed $p 128 > $out/$out.microsatellite.$p.128.wig
  # Let's make this at 100% scale
  java -jar -Xmx1g $tools/T2T-Polish/paf_util/bedToFixedWig.jar $out/$out.microsatellite.$p.128.bed $p 128 5 > $out/$out.microsatellite.$p.128.wig

  echo "Filter regions > 80%"
  awk -v n=80 '$NF>n' $out/$out.microsatellite.$p.128.bed | bedtools merge -i - > $out/$out.microsatellite.$p.128.gt80.bed

done

module load ucsc

set -x

for p in GC AT GA TC
do
  echo "Convert wig to bw"
  wigToBigWig -clip $out/$out.microsatellite.$p.128.wig $fa.fai $out/$out.microsatellite.$p.128.bw
done
