#!/bin/bash

if [[ "$#" -lt 4 ]]; then
  echo "Usage: ./issues.sh in.paf name ver platform"
  echo "  in.paf: input paf file"
  echo "  name:   name to appear on the .wig files"
  echo "  ver:    assembly version"
  echo "  platform: HiFi or ONT"
  echo
  echo "Required: pattern/ver/, ver.bed, ver.error.bed, ver.exclude.bed, ver.telo.bed"
  exit -1
fi
in_paf=$1
name=$2
in=${in_paf/.paf/}
ver=$3
platform=$4

set -x
set -o pipefail

if [[ ! -s $ver.bed ]]; then
  echo "Linking files in ../"
  [[ -s ../$ver.bed ]] && ln -s ../$ver.bed
  [[ -s ../$ver.error.bed ]] && ln -s ../$ver.error.bed
  [[ -s ../$ver.exclude.bed ]] && ln -s ../$ver.exclude.bed
  [[ -s ../$ver.telo.bed ]] && ln -s ../$ver.telo.bed
  [[ -e ../pattern ]] && ln -s ../pattern
fi

# Collect wig files
$tools/T2T-Polish/coverage/collect_summary.sh $in_paf "$name"

# Collect coverage statistics
$tools/T2T-Polish/coverage/coverage_stat.sh $in.cov.wig

# Collect low supportive regions
$tools/T2T-Polish/coverage/low_support.sh $in_paf $ver $platform

