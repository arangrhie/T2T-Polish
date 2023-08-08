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

	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > index.jid
	index=`cat index.jid`
	echo $index
	extra="--dependency=$index"
fi

if [[ ! -s bam.list ]]; then
  cpus=24
  mem=120g
  name=map.$out
  script=$PIPELINE/bwa/bwa.sh
  args="$ref $fastq_map $out"
  partition=norm
  walltime=2-0
  path=`pwd`
  log=logs/$name.%A_%a.log

  LEN=`wc -l $fastq_map | awk '{print $1}'`
  extra="--gres=lscratch:700 $extra"
  if [[ -z $line_num ]]; then
    extra="$extra --array=1-$LEN"
  else
    extra="$extra --array=$line_num"
  fi

  echo "\
  sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
  sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > bwa.jid
  wait_for="--dependency=afterok:"`cat bwa.jid`

  for i in $(seq 1 $LEN)
  do
    echo "$out.$i.bam" >> bam.list
    echo "$out.$i.dedup.bam" >> bam.dedup.list
  done

else
  wait_for=""
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

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args

args="$out.dedup bam.dedup.list"
echo "
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' >> mrg.jid


