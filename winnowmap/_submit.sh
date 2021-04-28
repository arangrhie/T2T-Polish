#!/bin/bash

if [[ "$#" -lt 3 ]] ; then
  echo "Usage: ./_submit.sh ref.fasta prefix map [wm_opt]"
  echo "  Required: input.fofn"
  exit -1
fi

ref=$1
prefix=$2
map=$3
wm_opt=$4

PIPELINE=$tools/T2T-Polish/winnowmap

if ! [[ -e input.fofn ]]; then
  echo "No input.fofn found. Exit."
  exit -1
fi

path=`pwd`
mkdir -p logs

if ! [[ -s repetitive_k15.txt ]]; then
  cpus=12
  mem=24g
  name=init.$prefix
  script=$PIPELINE/init.sh
  args="$ref"
  partition=quick
  walltime=30:00
  log=logs/$name.%A.log

  echo "
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > init.jid

  jid=`cat init.jid`
  extra="--dependency=afterok:$jid"
fi

cpus=24
mem=60g
name=map.$prefix
script=$PIPELINE/map.sh
args="$ref $map $wm_opt"
partition=norm
walltime=2-0

LEN=`wc -l input.fofn | awk '{print $1}'`

extra="$extra --array=1-$LEN" # include the init.jid
log=logs/$name.%A_%a.log

echo "\
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > map.jid

jid=`cat map.jid`

cpus=24
mem=60g
name=merge.$prefix
script=$PIPELINE/merge.sh
args="$prefix"
partition=norm
walltime=1-0
path=`pwd`
extra="--dependency=afterok:$jid"
log=logs/$name.%A.log

echo "\
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > merge.jid

jid=`cat merge.jid`
extra="--dependency=afterok:$jid"

sh ~/codes/_submit_norm.sh 12 32g tdf.$prefix $tools/IGVTools/to_tdf.sh "$prefix.bam $ref" $extra

sh ~/codes/_submit_norm.sh 12 8g filt.$prefix $PIPELINE/filt.sh $prefix.bam $extra > filt.jid

jid=`tail -n1 filt.jid`
extra="--dependency=afterok:$jid"

sh ~/codes/_submit_norm.sh 12 32g tdf.$prefix.pri $tools/IGVTools/to_tdf.sh "$prefix.pri.bam $ref" $extra

sh ~/codes/_submit_norm.sh 12 8g sam2paf.$prefix $PIPELINE/../coverage/sam2paf.sh "$prefix.pri.bam $prefix.pri.paf" $extra

