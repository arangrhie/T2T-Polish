#!/bin/bash

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./step1.sh ref.fa in.bam SAMPLE"
  echo "  ref.fa  : reference sequence"
  echo "  in.bam  : read alignment"
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

module load deepvariant/1.5.0 || exit 1
module load parallel

wd=$PWD

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID}
echo "Working temporarily in $PWD"

# assign REF to actual path of file, allowing for relative paths from the original working directory:
if [[ "$1" != "${1#/}" ]]; then
   REF=`realpath $1`
else
   REF=`realpath $wd/$1`
fi

if [[ "$2" != "${2#/}" ]]; then
   BAM=`realpath $2`
else
   BAM=`realpath $wd/$2`
fi
SAMPLE=$3
MODE=`cat $wd/MODE`
OUT=dv_$MODE/examples # written in /lscratch by default
N_SHARD=`cat $wd/N_SHARD`

mkdir -p $OUT logs

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
echo "
make_examples \\
  --mode calling   \\
  --ref "${REF}"   \\
  --reads "${BAM}" \\ 
  --examples $OUT/tfrecord@${N_SHARD}.gz $extra_args \\
  --sample_name "${SAMPLE}" \\
  --task {}
"

seq 0 $((N_SHARD-1)) \
  | parallel -j ${SLURM_CPUS_PER_TASK} --eta --halt 2 \
  --joblog "logs/log" --res "logs" \
  make_examples    \
    --mode calling \
    --ref "${REF}" \
    --reads "${BAM}" \
    --examples $OUT/tfrecord@${N_SHARD}.gz $extra_args \
    --sample_name "${SAMPLE}" \
    --task {} \
    || exit 1

mkdir -p $wd/$OUT
cp -r $OUT/* $wd/$OUT/
cp -r logs $wd/dv_$MODE/
