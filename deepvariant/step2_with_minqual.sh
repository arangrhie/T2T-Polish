#! /bin/bash

echo "Usage: ./step2_with_minqual.sh <minimum quality>"

MINQUAL=$1

if [[ ! -s MODE_MQ$MINQUAL ]]; then
  echo "Requires MODE_MQ[MINQUAL] text file to be present. Exit."
  exit -1
fi

module load deepvariant/1.5.0 || exit 1

set -o pipefail
set -e

MODE=`cat MODE_MQ$MINQUAL`
OUT=dv_$MODE\_MQ$MINQUAL
mode=`echo $MODE | awk '{print tolower($0)}'`
MODEL="/opt/models/$mode/model.ckpt"
N_SHARD=`cat N_SHARD_MQ$MINQUAL` # Must match as in step1
CALL_VARIANTS_OUTPUT="$OUT/call_variants_output.tfrecord.gz"

set -x

call_variants \
  --outfile "${CALL_VARIANTS_OUTPUT}" \
  --examples "$OUT/examples/tfrecord@${N_SHARD}.gz" \
  --checkpoint "${MODEL}" \
  || exit 1
