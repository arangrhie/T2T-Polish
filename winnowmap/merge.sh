#!/bin/bash

if [ -z $1 ]; then
  echo "Usage: merge.sh <out_prefix>"
  exit -1
fi

prefix=$1

module load samtools # load v1.15.1 or higher
echo "\
samtools merge -O bam -@$SLURM_CPUS_PER_TASK $prefix.bam *.sort.bam"
samtools merge -O bam -@$SLURM_CPUS_PER_TASK $prefix.bam *.sort.bam
echo "\
samtools index $prefix.bam"
samtools index $prefix.bam
