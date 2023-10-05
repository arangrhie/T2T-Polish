#! /bin/bash

echo "./step3_with_minqual.sh minqual"

set -o pipefail
set -e

module load deepvariant/1.5.0

export MINQUAL=$1

REF=`cat REF_MQ$MINQUAL`
MODE=`cat MODE_MQ$MINQUAL`
OUT=dv_$MODE\_MQ$MINQUAL
CALL_VARIANTS_OUTPUT="$OUT/call_variants_output.tfrecord.gz"

set -x
postprocess_variants \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "$OUT/$OUT.vcf.gz"

