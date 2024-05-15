#! /bin/bash

if [[ "$#" -lt 3 ]]; then
  echo "Usage: \$tools/T2T-Polish/deepvariant/_submit_ont_r9_pepper_margin_dv.sh bam ref mq [wait-for-jid]"
  echo "bam           input bam"
  echo "ref           reference fasta. REQUIRES .fai file in the same place"
  echo "mq            minimum mapping quality requirements. Use 0 for all-to-dip, 1 for all-to-hap alignment."
  echo "wait-for-jid  OPTIONAL. slurm job id to wait for."
  exit -1
fi

bam=$1
ref=$2
mq=$3
wait_for=$4

NR=`wc -l $ref.fai | awk '{print $1}'`

cpus=24
mem=48g
name=ont_dv.$bam
script=$tools/T2T-Polish/deepvariant/ont_r9_pepper_margin_dv.sh
args="$bam $ref $mq"
partition=gpu
walltime=2:00:00
path=`pwd`
extra="--array=1-$NR"

if ! [[ "$wait_for" == "" ]]; then
  extra="$extra --dependency=afterok:$wait_for"
fi

mkdir -p logs
log=logs/$name.%A_%a.log

set -x
sbatch -J $name \
  --cpus-per-task=$cpus --mem=$mem --partition=$partition \
  -D $path \
  --gres=gpu:k80:4,lscratch:10 \
  $extra --time=$walltime \
  --error=$log \
  --output=$log \
  $script $args > ont_dv.jid
set +x
cat ont_dv.jid

## Merge
wait_for=`cat ont_dv.jid | tail -n1`

mem=8g
name=ont_dv_mrg.$bam
script=$tools/T2T-Polish/deepvariant/merge_per_chr_vcfs.sh
args="dv_ONT_R9_ dv_ONT_R9"
if [[ $mq -eq -1 ]]; then
  args="${args}"
else
  args="${args}_MQ$mq"
fi
partition=quick
walltime=30:00
extra="--dependency=afterany:$wait_for" # afterany because dv reports "FAILED" for empty VCFs
log=logs/$name.%A.log

set -x
sbatch -J $name \
  --cpus-per-task=$cpus --mem=$mem --partition=$partition \
  -D $path \
  $extra --time=$walltime \
  --error=$log \
  --output=$log \
  $script $args > ont_dv.mrg.jid
set +x
cat ont_dv.mrg.jid

