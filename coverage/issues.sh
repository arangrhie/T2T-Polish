#!/bin/bash

if [[ "$#" -lt 5 ]]; then
  echo "Usage: ./issues.sh in.paf name ver platform pattern"
  echo "  in.paf   input paf file"
  echo "  name     name to appear on the .wig files"
  echo "  ver      assembly version"
  echo "  platform HiFi, ONT or Hybrid"
  echo "  pattern  path to pattern folder"
  echo
  echo "Required: pattern/ver/, pattern/ver.bed, pattern/ver.error.bed, pattern/ver.exclude.bed, pattern/ver.telo.bed"
  exit -1
fi
in_paf=$1
name=$2
ver=$3
platform=$4
pattern=$5

ln $in_paf
in_paf=`basename $in_paf`
in=${in_paf/.paf/}

if [[ ! -s $in_paf ]]; then
  echo "No $in_paf found. Exit."
  exit -1
fi

if [[ ! -d $pattern/$ver || ! -s $pattern/$ver.exclude.bed || ! -s $pattern/$ver.telo.bed || ! -s $pattern/$ver.bed ]]; then
  echo "No $pattern/$ver found. Exit. Run T2T-Polish/coverage/init.sh under pattern folder first."
  exit -1
fi

if [[ ! -s $pattern/$ver.error.bed ]]; then
  echo "No $pattern/$ver.error.bed found. Exit. Link the merged _only.bed file from Merqury as $pattern/$ver.error.bed."
  exit -1
fi

ln -sf $pattern pattern

set -x
set -o pipefail

# Collect wig files
$tools/T2T-Polish/coverage/collect_summary.sh $in_paf "$name"

# Collect coverage statistics
$tools/T2T-Polish/coverage/coverage_stat.sh $in.cov.wig pattern/$ver.exclude.bed

# Collect low supportive regions
$tools/T2T-Polish/coverage/low_support.sh $in_paf $ver $platform $pattern

awk -v OFS='\t' -v FS='\t' 'NR==FNR { mapping[$1] = $3; next } \
  {if ($3 > mapping[$1]) { $3=mapping[$1]; $8=mapping[$1];} print}' \
  $pattern/$ver.bed $in.issues.bed > $in.issues.fm.bed

module load ucsc
wigToBigWig -clip $in.cov.wig pattern/$ver.fa.fai $in.cov.bw
bedToBigBed -type=bed9 $in.issues.fm.bed pattern/$ver.fa.fai $in.issues.bb

