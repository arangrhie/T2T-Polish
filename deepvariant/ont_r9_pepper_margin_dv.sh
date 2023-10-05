#!/bin/sh

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./ont_r9_pepper_margin_dv.sh bam ref mq [region]"
  echo "  bam     input bam"
  echo "  ref     reference fasta. REQUIRES .fai file in the same place"
  echo "  mq      minimum mapping quality requirements. Use 0 for all-to-dip, 1 for all-to-hap alignment."
  echo "  region  OPTIONAL. Looks for SLURM_ARRAY_JOB_ID first, if not found, tries to use this."
  exit -1
fi

BAM=$1
REF=$2
MQ=$3

THREADS=$SLURM_CPUS_PER_TASK

i=$SLURM_ARRAY_JOB_ID
REGION=`sed -n ${i}p $REF.fai | cut -f1`
if [[ "$REGION" == "" ]]; then
  REGION=$4
fi

OUTPUT_DIR=dv_ONT_R9_MQ${MQ}_$REGION

ulimit -u 8192

set -e
set -o pipefail

module load pepper_deepvariant/0.8

set -x

run_pepper_margin_deepvariant call_variant \
  -b $BAM -f $REF -o $OUTPUT_DIR \
  --pepper_min_mapq $MQ --dv_min_mapping_quality $MQ \
  -t $THREADS -r $REGION --ont_r9_guppy5_sup --gpu

