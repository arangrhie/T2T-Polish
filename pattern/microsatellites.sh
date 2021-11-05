#!/bin/bash

if [[ -z $1 ]]; then
  echo "Usage: ./microsatellites.sh in.fasta"
  exit -1
fi
fa=$1

for p in ga gc at
do
  echo "Collect $p"
  $tools/seqrequester/build/bin/seqrequester microsatellite -prefix microsatellite -window 128 -$p $fa
done

for p in GC AT GA TC
do
  echo "Convert to bed: $p"
  java -jar -Xmx1g $tools/T2T-Polish/paf_util/bedToFixedWig.jar microsatellite.$p.128.bed $p 128 > microsatellite.$p.128.wig

  echo "Filter regions > 80%"
  awk -v n=80 '$NF>n' microsatellite.$p.128.bed > microsatellite.$p.128.gt80.bed
done
