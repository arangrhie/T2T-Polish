#!/bin/bash

if [[ "$#" -lt 3 ]] ; then
  echo "Usage: ./_submit.sh ref.fasta prefix map [wm_opt]"
  echo "  ref.fasta  reference to align"
  echo "  prefix     output prefix"
  echo "  map        mapping mode. map-pb for HiFi, map-ont for ONT, map-pb-clr for CLR."
  echo "  wm_opt     winnowmap optional arguments. i.e. -y for adding methylation tags."
  echo "  Required: input.fofn"
  exit -1
fi

ref=$1
prefix=$2
map=$3
wm_opt=$4
extra=$5 # "--dependency=afterok:$jid"

PIPELINE=$tools/T2T-Polish/winnowmap

echo $wn_opt


if ! [[ -e input.fofn ]]; then
  echo "No input.fofn found. Exit."
  exit -1
fi

path=`pwd`
mkdir -p logs

ln -s $ref
ref=`basename $ref`

set -o pipefail


if ! [[ -s repetitive_k15.txt ]]; then
  cpus=12
  mem=24g
  partition=quick
  walltime=30:00
  name=init.$prefix
  log=logs/$name.%A.log
  script=$PIPELINE/init.sh
  args="$ref"

  set -x
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D `pwd` $extra --time=$walltime --error=$log --output=$log $script $args > init.jid
  set +x

  jid=`cat init.jid`
  echo $jid

  extra="--dependency=afterok:$jid"
fi

cpus=24
mem=120g
partition=norm
walltime=2-0
name=map.$prefix
log=logs/$name.%A_%a.log
script=$PIPELINE/map.sh
# args="$ref $map $wm_opt"

echo $args
LEN=`wc -l input.fofn | awk '{print $1}'`
arr=""
for i in $(seq 1 $LEN)
do
  reads=`sed -n ${i}p input.fofn`

  out=`basename $reads`
  out=`echo $out | sed 's/.gz$//g'`
  out=`echo $out | sed 's/.fasta$//g' | sed 's/.fa$//g'`
  out=`echo $out | sed 's/.fastq$//g' | sed 's/.fq$//g'`
  out=$out.$i

  if ! [[ -s $out.sort.bam ]] ; then
    arr="${arr}${i}_"
  fi
done

# bash does not understand commas (,), so let's use _
if [[ "$arr" == "" ]]; then
  echo "Found all *.sort.bam. Skip mapping"
else
  # add 900g local sractch, include job dependency to init.jid
  arr=`echo $arr | sed 's/_/,/g'`
  extra="$extra --gres=lscratch:900 --array=$arr" 

  set -x
  # sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D `pwd` $extra --time=$walltime --error=$log --output=$log $script $args > map.jid
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D `pwd` $extra --time=$walltime --error=$log --output=$log $script $ref $map "$wm_opt" > map.jid
  set +x
  cat map.jid
fi

# Merge
cpus=48
mem=60g
partition=norm
walltime=1-0
name=merge.$prefix
script=$PIPELINE/merge.sh
args="$prefix"

jid=`cat map.jid`
extra="--dependency=afterok:$jid"
log=logs/$name.%A.log

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D `pwd` $extra --time=$walltime --error=$log --output=$log $script $args > merge.jid
set +x

cpus=12
mem=8g
name=filt.$prefix
log=logs/$name.%A.log
script=$PIPELINE/filt.sh
args="$prefix.bam"

jid=`cat merge.jid`
extra="--dependency=afterok:$jid"
set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D `pwd` $extra --time=$walltime --error=$log --output=$log $script $args > filt.jid
set +x
cat filt.jid

name=sam2paf.$prefix
log=logs/$name.%A.log
script=$PIPELINE/../coverage/sam2paf.sh
args="$prefix.pri.bam $prefix.pri.paf"

jid=`tail -n1 filt.jid`
extra="--dependency=afterok:$jid"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem --partition=$partition -D `pwd` $extra --time=$walltime --error=$log --output=$log $script $args
set +x
