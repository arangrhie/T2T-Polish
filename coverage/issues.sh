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
in=${in_paf/.paf/}
ver=$3
platform=$4
pattern=$5

ln -sf $pattern pattern

set -x
set -o pipefail

# Collect wig files
$tools/T2T-Polish/coverage/collect_summary.sh $in_paf "$name"

# Collect coverage statistics
$tools/T2T-Polish/coverage/coverage_stat.sh $in.cov.wig pattern/$ver.exclude.bed

# Collect low supportive regions
$tools/T2T-Polish/coverage/low_support.sh $in_paf $ver $platform $pattern

