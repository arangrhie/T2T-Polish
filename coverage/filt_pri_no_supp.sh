#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./filt_pri_no_supp.sh in.bam"
  echo "  output: in.pri.bam and in.pri.paf"
  exit -1
fi

BAM=$1
PRI=${BAM/.bam/.pri.bam}
PAF=${PRI/.bam/.paf}

CPU=$SLURM_CPUS_PER_TASK

samtools view -F0x900 -hb -@$CPU $BAM > $PRI
samtools index $PRI
$tools/T2T-Polish/coverage/sam2paf.sh $PRI $PAF

