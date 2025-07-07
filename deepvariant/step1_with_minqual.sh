#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./step1_with_minqual.sh SAMPLE MINQUAL"
  echo "  SAMPLE  : SAMPLE name to be present in the final VCF"
  echo "  MINQUAL  : Value of --min-mapping-quality to pass to make_examples and use in output directory name"
  echo 
  echo "Requires REF_MQ[MINQUAL] BAM_MQ[MINQUAL] MODE_MQ[MINQUAL] and N_SHARD_MQ[MINQUAL] text files."
  echo "  MODE_MQ[MINQUAL]    : WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "                        MODE specific paramters are set accordingly as in https://github.com/google/deepvariant/blob/r1.5/scripts/run_deepvariant.py"
  echo "                        For WES, settings are the same. Use WGS instead."
  echo "  N_SHARD_MQ[MINQUAL] : Num. of threads / shards to use"
  exit -1
fi

set -e
set -o pipefail

# Update to DeepVariant v1.6.1 on Mar. 10 2025
module load deepvariant/1.6.1 || exit 1
module load parallel

set -x

wd=$PWD
SAMPLE=$1
MINQUAL=$2

if [[ $MINQUAL -eq -1 ]]; then
  echo "Use pre-set options in DV for --min_mapping_quality."
  minqual=""
else
  echo "Use --min_mapping_quality $MINQUAL"
  minqual="--min_mapping_quality $MINQUAL"
fi

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID}
echo "Working temporarily in $PWD"

REF=`cat ${wd}/REF_MQ${MINQUAL}`
BAM=`cat ${wd}/BAM_MQ${MINQUAL}`
MODE=`cat ${wd}/MODE_MQ${MINQUAL}`
N_SHARD=`cat ${wd}/N_SHARD_MQ${MINQUAL}`
OUT=dv_${MODE}_MQ${MINQUAL}/examples # written in /lscratch by default

mkdir -p $OUT logs ${wd}/logs-parallel-$SLURM_JOB_ID

extra_args=""

if [[ $MODE == "WGS" ]]; then
  extra_args="--channels insert_size $minqual"
elif [[ $MODE == "PACBIO" ]]; then
  extra_args="--add_hp_channel --alt_aligned_pileup diff_channels --max_reads_per_partition 600 $minqual --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 199 --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12"
elif [[ $MODE == "ONT_R104" ]]; then
  extra_args="--add_hp_channel --alt_aligned_pileup diff_channels --max_reads_per_partition 600 $minqual --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 199 --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12 --vsc_min_fraction_snps 0.08"
elif [[ $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  extra_args="$minqual"
else
  echo "Unknown $MODE provided. Exit."
  exit -1
fi

echo "make_examples with parallel in $MODE mode"

# GVCF_TFRECORDS="${OUT}/examples.gvcf.tfrecord@${N_SHARDS}.gz"
# VCF_TFRECORDS="${OUT}/examples.tfrecord@${N_SHARDS}.gz"

seq 0 $(($N_SHARD-1)) \
  | parallel -j ${SLURM_CPUS_PER_TASK} --eta --halt 2 \
  --joblog "$wd/logs-parallel-$SLURM_JOB_ID/log" --res "$wd/logs-parallel-$SLURM_JOB_ID" \
    make_examples      \
      --mode calling   \
      --ref ${REF}   \
      --reads ${BAM} \
      --examples ${OUT}/examples.tfrecord@${N_SHARD}.gz \
      --gvcf ${OUT}/examples.gvcf.tfrecord@${N_SHARD}.gz    \
      --sample_name $SAMPLE     \
      $extra_args \
      --task {}

mkdir -p ${wd}/$OUT &&
mv ${OUT} ${wd}/dv_${MODE}_MQ${MINQUAL}/ &&
touch ${wd}/deepvariant.step1.done || exit 1