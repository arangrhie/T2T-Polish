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
  partition=quick
  walltime=30:00
  name=init.$prefix
  log=logs/$name.%A.log
  script=$PIPELINE/init.sh
  args="$ref"

  echo "
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > init.jid

  jid=`cat init.jid`
  extra="--dependency=afterok:$jid"
fi

cpus=24
mem=60g
partition=norm
walltime=2-0
name=map.$prefix
log=logs/$name.%A_%a.log
script=$PIPELINE/map.sh
args="$ref $map $wm_opt"

LEN=`wc -l input.fofn | awk '{print $1}'`
extra="$extra --array=1-$LEN" # include job dependency to init.jid

echo "\
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > map.jid

cpus=24
mem=60g
partition=norm
walltime=1-0
name=merge.$prefix
script=$PIPELINE/merge.sh
args="$prefix"

jid=`cat map.jid`
extra="--dependency=afterok:$jid"
log=logs/$name.%A.log

echo "\
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > merge.jid

cpus=12
mem=8g
name=filt.$prefix
logs=logs/$name.%A.log
script=$PIPELINE/filt.sh
args="$prefix.bam"

jid=`cat merge.jid`
extra="--dependency=afterok:$jid"
echo "\
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args > filt.jid

name=sam2paf.$prefix
logs=logs/$name.%A.log
script=$PIPELINE/../coverage/sam2paf.sh
args="$prefix.pri.bam $prefix.pri.paf"

jid=`tail -n1 filt.jid`
extra="--dependency=afterok:$jid"

echo "\
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D $path $extra --time=$walltime --error=$log --output=$log $script $args

