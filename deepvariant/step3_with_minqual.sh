#! /bin/bash

echo "./step3_with_minqual.sh minqual"

set -o pipefail
set -e

module load deepvariant/1.6.1 || exit 1

export MINQUAL=$1

REF=`cat REF_MQ${MINQUAL}`
MODE=`cat MODE_MQ${MINQUAL}`
OUT=dv_${MODE}_MQ${MINQUAL}
N_SHARD=`cat N_SHARD_MQ${MINQUAL}`
CALL_VARIANTS_OUTPUT="${OUT}/call_variants_output.tfrecord.gz"
# GVCF_TFRECORDS="${OUT}/examples/examples.gvcf.tfrecord@${N_SHARD}.gz"
GVCF_TFRECORDS=${OUT}/examples/examples.gvcf.tfrecord@${N_SHARD}.gz
CPU=12
SAMPLE=`cat SAMPLE_MQ${MINQUAL}`

set -x
mkdir -p $OUT/nonvariant_site_tfrecord_path/

postprocess_variants \
  --ref ${REF} \
  --infile ${CALL_VARIANTS_OUTPUT} \
  --outfile ${OUT}.${SAMPLE}.vcf.gz \
  --gvcf_outfile ${OUT}.${SAMPLE}.gvcf.gz \
  --nonvariant_site_tfrecord_path $GVCF_TFRECORDS \
  --cpus $CPU \
  --sample_name ${SAMPLE} &&
touch deepvariant.step3.done || exit -1
set +x

echo "deepvariant done"


