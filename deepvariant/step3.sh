#! /bin/bash

echo "./step3.sh"

set -o pipefail
set -e

module load deepvariant/1.5.0

REF=`cat REF`
MODE=`cat MODE`
OUT=dv_$MODE
CALL_VARIANTS_OUTPUT="$OUT/call_variants_output.tfrecord.gz"

set -x
postprocess_variants \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "$OUT/$OUT.vcf.gz"

