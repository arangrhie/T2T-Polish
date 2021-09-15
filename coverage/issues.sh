#!/bin/bash

if [[ "$#" -lt 4 ]]; then
  echo "Usage: ./issues.sh in.paf name ver platform"
  echo "  in.paf: input paf file"
  echo "  name:   name to appear on the .wig files"
  echo "  ver:    assembly version"
  echo "  platform: HiFi or ONT"
  exit -1
fi
in_paf=$1
name=$2
in=${in_paf/.paf/}
ver=$3
platform=$4

# Collect wig files
$tools/T2T-Polish/coverage/collect_summary.sh $in_paf $name

# Collect coverage statistics
$tools/T2T-Polish/coverage/coverage_stat.sh $in.cov.wig

# Collect low supportive regions
$tools/T2T-Polish/coverage/low_support.sh $in_paf $ver $platform

