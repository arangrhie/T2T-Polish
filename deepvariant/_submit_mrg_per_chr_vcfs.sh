#!/bin/bash

echo "Usage: ./_submit_merge_per_chr_vcfs.sh in_dir out [jid]"
echo "  in_dir  input dir prefix. ex. dv_ONT_R9_"
echo "  out     output prefix. out.vcf.gz and out.visual_report.html will be generated."
echo "  jid     job id to wait for."

if [[ "$#" -lt 2 ]]; then
  exit -1
fi

in_dir=$1
out=$2
wait_for=$3
extra=""

cpus=12
mem=8g
name=dv.mrge
script=$tools/T2T-Polish/deepvariant/merge_per_chr_vcfs.sh
args="$in_dir $out"
if ! [[ -z $wait_for ]]; then
  extra="--dependency=afterok:$wait_for"
fi

mkdir -p logs

log=logs/$name.%A.log

set -x
sbatch --partition=quick \
      -D `pwd` \
      --job-name=$name --time=240 \
      $extra \
      --cpus-per-task=$cpus --mem=$mem \
      --error=$log \
      --output=$log \
      $script $args

