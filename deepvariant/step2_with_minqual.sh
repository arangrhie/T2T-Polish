#! /bin/bash

echo "Usage: ./step2_with_minqual.sh <minimum quality>"

MINQUAL=$1

if [[ ! -s MODE_MQ$MINQUAL ]]; then
  echo "Requires MODE_MQ[MINQUAL] text file to be present. Exit."
  exit -1
fi

# Update to DeepVariant v1.6.1 on Sep. 26 2024
module load deepvariant/1.6.1 || exit 1

ulimit -u 32768
set -o pipefail
set -e

MODE=`cat MODE_MQ$MINQUAL`
mode=`echo $MODE | awk '{print tolower($0)}'`
OUT=dv_${MODE}_MQ${MINQUAL}/examples
mkdir -p $OUT
MODEL="/opt/models/$mode"
N_SHARD=`cat N_SHARD_MQ$MINQUAL` # Must match as in step1
# CALL_VARIANTS_OUTPUT="dv_$MODE\_MQ$MINQUAL/call_variants_OUTput.tfrecord.gz"
CALL_VARIANTS_OUTPUT="$OUT/../call_variants_output.tfrecord.gz"

if ! ls $OUT/../*.gz 1> /dev/null 2>&1; then
  set -x

  call_variants \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --writer_threads=4 \
    --examples "$OUT/examples.tfrecord@${N_SHARD}.gz" \
    --checkpoint "${MODEL}"

  set +x
fi

# Lots of silently dying threads (not detectable with set -e, exit code is 0) found in the log.
# Let's check the actual file size of the OUTput.
# We see 0 file size matching the # threads died in the log.
ls -l $OUT/../call_variants_output*.tfrecord.gz

incompleteFiles=0
for cvo in $(ls $OUT/../call_variants_output*.tfrecord.gz)
do
  if [[ ! -s $cvo ]]; then
    echo "0 size file $cvo detected."
    incompleteFiles=$((incompleteFiles + 1))
  fi
done

echo "Incomplete Files: $incompleteFiles"
if [[ $incompleteFiles -eq 0 ]]; then
  touch deepvariant.step2.done
else
  echo "Try one more time in 5 sec..."
  sleep 5
  set -x
  call_variants \
  --outfile "${CALL_VARIANTS_OUTPUT}" \
  --writer_threads=4 \
  --examples $OUT/examples.tfrecord@${N_SHARD}.gz \
  --checkpoint "${MODEL}"
  set +x
  
  ls -l dv_${MODE}_MQ${MINQUAL}/call_variants_output*.tfrecord.gz

  incompleteFiles=0
  for cvo in $(ls $OUT/../call_variants_output*.tfrecord.gz)
  do
    if [[ ! -s $cvo ]]; then
      echo "0 size file $cvo detected."
      incompleteFiles=$((incompleteFiles + 1))
    fi
  done

  echo "Incomplete Files: $incompleteFiles"
  if [[ $incompleteFiles -eq 0 ]]; then
    touch deepvariant.step2.done
  else
    "Error still persists. Giving up."
    exit -1
  fi
fi