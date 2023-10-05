#! /bin/bash

echo "./step3.sh CHROM"

set -o pipefail
set -e

module load deepvariant/1.5.0-deeptrio
module load parallel

CHROM=$1
REF=`cat REF`
MODE=`cat MODE`
PARENT1SAMPLE=`cat PARENT1SAMPLE`
PARENT2SAMPLE=`cat PARENT2SAMPLE`
CHILDSAMPLE=`cat CHILDSAMPLE`
N_SHARDS=`cat N_SHARD`
OUT=$CHROM\_$MODE\_output
CALL_VARIANTS_OUTPUT="$OUT/call_variants_output.tfrecord.gz"

declare -A samples=( [parent1]=$PARENT1SAMPLE [parent2]=$PARENT2SAMPLE [child]=$CHILDSAMPLE )

for n in parent1 parent2 child ; do
    postprocess_variants \
        --ref $REF \
        --infile $CHROM\_WGS_output/call_variants_output_${n}.tfrecord.gz \
        --outfile $CHROM\_WGS_output/${samples[$n]}.output.vcf.gz \
        --nonvariant_site_tfrecord_path $CHROM\_WGS_output/nvgvcf/gvcf_${n}.tfrecord@${N_SHARDS}.gz \
        --gvcf_outfile $CHROM\_WGS_output/${samples[$n]}.g.vcf.gz
done

