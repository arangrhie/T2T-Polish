#! /bin/sh

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./_submit_bwa.sh <ref.fasta> <fastq.map> <out> [line_num]"
  exit -1
fi

ref=$1
fastq_map=$2
out=$3
line_num=$4

mkdir -p logs
PIPELINE=$tools/T2T-Polish

set -e
set -o pipefail

if [[ ! -s $ref.bwt ]]; then
	cpus=4
	mem=10g
	name=index.$out
	script=$PIPELINE/bwa/bwa_index.sh
	args="$ref"
	partition=quick
	walltime=4:00:00
	path=`pwd`
	log=logs/$name.%A.log

  set -x
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > index.jid
  set +x
	index=`cat index.jid`
	echo "Indexing: $index"
	extra="--dependency=$index"
fi

LEN=`wc -l $fastq_map | awk '{print $1}'`
arr=""
for i in $(seq 1 $LEN)
do
  if [[ ! -s $out.$i.dedup.bam ]]; then
    arr="${arr}${i}_"
  fi
done

if [[ -v "$arr" ]]; then
  echo "Found all *.dedup.bam. Skip mapping"
  extra=""
  wait_for=""
else
  arr=`echo $arr | sed 's/_/,/g'`
  cpus=24
  mem=120g
  name=map.$out
  script=$PIPELINE/bwa/bwa.sh
  args="$ref $fastq_map $out"
  partition=norm
  walltime=2-0
  path=`pwd`
  log=logs/$name.%A_%a.log

  extra="--gres=lscratch:900 $extra"
  if [[ -z $line_num ]]; then
    extra="$extra --array=$arr"
  else
    extra="$extra --array=$line_num"
  fi

  # map
  set -x
  sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > bwa.jid
  set +x
  jid=`cat bwa.jid`
  echo "Map: $jid"
  wait_for="--dependency=afterok:$jid"

  for i in $(seq 1 $LEN)
  do
    echo "$out.$i.bam" >> bam.list
    echo "$out.$i.dedup.bam" >> bam.dedup.list
  done
fi

cpus=24
mem=48g
name=merge.$out
script=$PIPELINE/bwa/merge.sh
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A.log
extra=$wait_for

args="$out bam.list"

set -x
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args
set +x

name=mrg.$out.dedup
args="$out.dedup bam.dedup.list"
set -x
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > mrg.jid
set +x
jid=`cat mrg.jid`
echo "Merge: $jid"
