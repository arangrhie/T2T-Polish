#!/bin/bash

if [[ -z $1 ]]; then
  echo "Usage: ./filt.sh in.bam"
  echo "  output: in.pri.bam"
  echo "          this bam file is filtered with -F0x104"
  exit -1
fi

in=$1
out=${in/.bam/.pri.bam}

module load samtools # load v1.15.1 or higher
echo "
samtools view -F0x104 -@$SLURM_CPUS_PER_TASK -hb $in > $out"
samtools view -F0x104 -@$SLURM_CPUS_PER_TASK -hb $in > $out

echo "
samtools index $out"
samtools index $out

