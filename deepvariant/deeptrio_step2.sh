#! /bin/bash

echo "Usage: ./deeptrio_step2.sh CHROM"
echo "  CHROM : reference chromsome for this calling step"

if [[ ! -s MODE ]]; then
  echo "Requires MODE text file to be present. Exit."
  exit -1
fi

module load deepvariant/1.5.0-deeptrio || exit 1
module load parallel || exit 1

set -o pipefail
set -e

CHROM=$1
MODE=`cat MODE`
OUT=$CHROM\_$MODE\_output
mode=`echo $MODE | awk '{print tolower($0)}'`
N_SHARDS=`cat N_SHARD` # Must match as in step1

set -x

for n in parent1 parent2 child ; do
    modeldir=`echo $n | sed 's/[12]//'`
    call_variants \
        --outfile $OUT/call_variants_output_${n}.tfrecord.gz \
        --examples $OUT/examples/make_examples_${n}.tfrecord@${N_SHARDS}.gz \
        --checkpoint /opt/models/deeptrio/wgs/$modeldir/model.ckpt
done

