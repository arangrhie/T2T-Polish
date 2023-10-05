#!/bin/bash

if [[ "$#" -lt 8 ]]; then
  echo "Usage: ./step1.sh ref.fa child.bam parent1.bam parent2.bam CHILDSAMPLE PARENT1SAMPLE PARENT2SAMPLE"
  echo "  ref.fa       : reference sequence"
  echo "  child.bam    : childs read alignment"
  echo "  parent1.bam  : parent1s read alignment"
  echo "  parent2.bam  : parent2s read alignment"
  echo "  CHILDSAMPLE  : child's sample name to be present in the final VCF"
  echo "  PARENT1SAMPLE  : parent1's sample name to be present in the final VCF"
  echo "  PARENT2SAMPLE  : parent2's sample name to be present in the final VCF"
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

module load deepvariant/1.5.0-deeptrio || exit 1
module load parallel
module load samtools

wd=$PWD

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID}
echo "Working temporarily in $PWD"

REF=$1
CHILDBAM=$2
PARENT1BAM=$3
PARENT2BAM=$4
CHILDSAMPLE=$5
PARENT1SAMPLE=$6
PARENT2SAMPLE=$7
CHROM=$8
MODE=`cat $wd/MODE`
OUT=$CHROM\_$MODE\_output
EXOUT=$OUT/examples # written in /lscratch by default
NVOUT=$OUT/nvgvcf
N_SHARDS=`cat $wd/N_SHARD`

mkdir -p $OUT $EXOUT $NVOUT

# extract chromosome alignments:
mkdir input
CHILDCHROMBAM=`echo $CHILDBAM | sed 's:.*/::'` 
samtools view -h $CHILDBAM $CHROM | samtools view -b > input/$CHILDCHROMBAM
samtools index input/$CHILDCHROMBAM
PARENT1CHROMBAM=`echo $PARENT1BAM | sed 's:.*/::'` 
samtools view -h $PARENT1BAM $CHROM | samtools view -b > input/$PARENT1CHROMBAM
samtools index input/$PARENT1CHROMBAM
PARENT2CHROMBAM=`echo $PARENT2BAM | sed 's:.*/::'` 
samtools view -h $PARENT2BAM $CHROM | samtools view -b > input/$PARENT2CHROMBAM
samtools index input/$PARENT2CHROMBAM

extra_args=""

if [[ $MODE == "WGS" ]]; then
  extra_args="--channels insert_size"
else
  echo "Haven\'t implemented mode $MODE. Exit."
  exit -1
fi

echo "make_examples with parallel in $MODE mode"
echo "
make_examples \\
  --mode calling   \\
  --ref "${REF}"   \\
  --reads_parent1 "input/${PARENT1CHROMBAM}" \\ 
  --sample_name_parent1 "${PARENT1SAMPLE}" \\ 
  --reads_parent2 "input/${PARENT2CHROMBAM}" \\ 
  --sample_name_parent2 "${PARENT2SAMPLE}" \\ 
  --reads "input/${CHILDCHROMBAM}" \\ 
  --sample_name "${CHILDSAMPLE}" \\
  --regions $CHROM \\
  --examples $EXOUT/make_examples.tfrecord@${N_SHARDS}.gz $extra_args \\
  --gvcf $NVOUT/gvcf.tfrecord@${N_SHARDS}.gz \\
  --pileup_image_height_child "60" \\
  --pileup_image_height_parent "40" \\
  --task {}
"

seq 0 $((N_SHARDS-1)) \
  | parallel -j ${SLURM_CPUS_PER_TASK} --eta --halt 2 \
  --joblog "$OUT/logs" --res "$OUT" \
  make_examples    \
    --mode calling \
    --ref "${REF}" \
    --reads_parent1 "input/${PARENT1CHROMBAM}" \
    --sample_name_parent1 "${PARENT1SAMPLE}" \
    --reads_parent2 "input/${PARENT2CHROMBAM}" \
    --sample_name_parent2 "${PARENT2SAMPLE}" \
    --reads "input/${CHILDCHROMBAM}" \
    --sample_name "${CHILDSAMPLE}" \
    --regions $CHROM \
    --examples $EXOUT/make_examples.tfrecord@${N_SHARDS}.gz $extra_args \
    --gvcf $NVOUT/gvcf.tfrecord@${N_SHARDS}.gz \
    --pileup_image_height_child "60" \
    --pileup_image_height_parent "40" \
    --task {}

mkdir -p $wd/$OUT
cp -r $OUT/* $wd/$OUT/
