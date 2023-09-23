#!/bin/sh

## merge_hybrid.sh

if [[ "$#" -lt 3 ]]; then
  echo "Usage: merge_hybrid.sh hifi.bam ilmn.bam out.bam"
  exit -1
fi

BAM_HIFI=$1
BAM_ILMN=$2
BAM_HYBR=$3

module load samtools
set -o pipefail
set -e
set -x

samtools merge -@$SLURM_CPUS_PER_TASK -O bam -o $BAM_HYBR $BAM_HIFI $BAM_ILMN
samtools index $BAM_HYBR


