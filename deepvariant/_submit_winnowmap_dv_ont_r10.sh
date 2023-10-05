#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: \$tools/T2T-Polish/deepvariant/_submit_winnowmap_dv.sh ref.fa out-prefix"
  echo "  ref         reference fasta. REQUIRES .fai file in the same place"
  echo "  out-prefix  output prefix"
  echo "  mq          minimum mapping quality requirements. Use 0 for all-to-dip, 1 for all-to-hap alignment."
  exit -1
fi

ref=$1
out=$2
mq=$3
map=map-ont

$tools/T2T-Polish/winnowmap/_submit.sh $ref $out $map
jid=`cat filt.jid | tail -n1`

$tools/T2T-Polish/deepvariant/_submit_deepvariant_with_minqual.sh $ref $out.pri.bam ONT_R10 ONT $mq $jid

