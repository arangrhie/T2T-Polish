#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./step1.sh SAMPLE"
  echo "  SAMPLE  : SAMPLE name to be present in the final VCF"
  echo 
  echo "Requires MODE and N_SHARD text files."
  echo "  MODE    : WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "            MODE specific paramters are set accordingly as in https://github.com/google/deepvariant/blob/r1.5/scripts/run_deepvariant.py"
  echo "            For WES, settings are the same. Use WGS instead."
  echo "  N_SHARD : Num. of threads / shards to use"
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

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID}
echo "Working temporarily in $PWD"

# Collect REF BAM MODE and N_SHARD
REF=`cat ${wd}/REF`
BAM=`cat ${wd}/BAM`
MODE=`cat ${wd}/MODE`
N_SHARD=`cat ${wd}/N_SHARD`
OUT=dv_$MODE/examples # written in /lscratch by default

mkdir -p $OUT logs $wd/logs-parallel-$SLURM_JOB_ID

extra_args=""

if [[ $MODE == "WGS" ]]; then
  extra_args="--channels insert_size"
elif [[ $MODE == "PACBIO" ]]; then
  extra_args="--add_hp_channel --alt_aligned_pileup diff_channels --max_reads_per_partition 600 --min_mapping_quality 1 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 199 --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12"
elif [[ $MODE == "ONT_R104" ]]; then
  extra_args="--add_hp_channel --alt_aligned_pileup diff_channels --max_reads_per_partition 600 --min_mapping_quality 5 --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 199 --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12 --vsc_min_fraction_snps 0.08"
elif [[ $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  extra_args=""
else
  echo "Unknown $MODE provided. Exit."
  exit -1
fi

echo "make_examples with parallel in $MODE mode"

GVCF_TFRECORDS="${OUT}/examples.gvcf.tfrecord@${N_SHARDS}.gz"
VCF_TFRECORDS="${OUT}/examples.tfrecord@${N_SHARDS}.gz"

seq 0 $((N_SHARD-1)) \
  | parallel -j ${SLURM_CPUS_PER_TASK} --eta --halt 2 \
  --joblog "$wd/logs-parallel-$SLURM_JOB_ID/log" --res "$wd/logs-parallel-$SLURM_JOB_ID" \
  make_examples      \
    --mode calling   \
    --ref   "${REF}" \
    --reads "${BAM}" \
    --examples $$VCF_TFRECORDS \
    --gvcf     $GVCF_TFRECORDS \
    --sample_name "${SAMPLE}"  \
    $extra_args \
    --task {} &&
mkdir -p $wd/dv_$MODE &&
cp -r $OUT $wd/dv_$MODE &&
touch $wd/deepvariant.step1.done || exit 1
