#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./map.sh ref map opt [i]"
  echo "Map reads with winnowmap."
  echo "Requires input.fofn and repetitive_k15.txt"
  echo "  ref    : ref.fa"
  echo "  map    : Provide one from map-ont map-pb map-pb-clr"
  echo "  opt    : additional options. add -y for adding methylation tags"
  echo "  [i]    : fq file to proceed in line i of input.fofn"
  echo
  echo "  output : input fq prefix appended with .sam and an indexed .sort.bam"
  exit -1
fi

cpus=$(($SLURM_CPUS_PER_TASK-2))
if [[ -z $cpus ]]; then
  cpus=12
fi

i=$SLURM_ARRAY_TASK_ID
if [[ -z $i ]]; then
  i=$4
fi

ref=$1
map=$2
opt=$3

tmp=/lscratch/${SLURM_JOB_ID}

reads=`sed -n ${i}p input.fofn`

out=`basename $reads`
out=`echo $out | sed 's/.gz$//g'`
out=`echo $out | sed 's/.fasta$//g' | sed 's/.fa$//g'`
out=`echo $out | sed 's/.fastq$//g' | sed 's/.fq$//g'`
out=$out.$i

set -e
set -o pipefail

if [[ -s $out.sort.bam ]] ; then
  echo "Found $out.sort.bam. Exit 0"
  exit 0
fi

module load winnowmap/2.03
module load samtools # load v1.15.1+

set -x
winnowmap --MD -W repetitive_k15.txt -ax $map -I12g $opt -t$cpus $ref $reads > $tmp/$out.sam

samtools sort -@$cpus -m2G -T $tmp/$out.tmp -O bam -o $tmp/$out.sort.bam $tmp/$out.sam

mv $tmp/$out.sort.bam ./

samtools index $out.sort.bam
set +x

