#!/bin/bash

if [[ "$#" -lt 3 ]]; then
  echo "Usage: \$tools/T2T-Polish/deepvariant/_submit_winnowmap_dv.sh ref.fa out-prefix mq [-y]"
  echo "  ref         reference fasta. REQUIRES .fai file in the same place"
  echo "  out-prefix  output prefix"
  echo "  mq          minimum mapping quality requirements. Use positive values for MQ filters. Use -1 to follow pre-set options in DV."
  echo "              e.g. 0 for all-to-dip, -1 for all-to-hap alignments."
  echo "  -y          OPTIONAL. set for passing methylation tags."
  exit -1
fi

ref=$1
out=$2
mq=$3
map="map-pb $4"

$tools/T2T-Polish/winnowmap/_submit.sh $ref $out "$map"
jid=`cat filt.jid | tail -n1`

if [[ "$mq" -eq -1 ]]; then
  $tools/T2T-Polish/deepvariant/_submit_deepvariant.sh $ref $out.pri.bam PACBIO HiFi $jid
else
  $tools/T2T-Polish/deepvariant/_submit_deepvariant_with_minqual.sh $ref $out.pri.bam PACBIO HiFi $mq $jid
fi

