#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: \$tools/T2T-Polish/deepvariant/_submit_winnowmap_dv_ont_r9.sh ref.fa out-prefix mq [jid]"
  echo "  ref         reference fasta. REQUIRES .fai file in the same place"
  echo "  out-prefix  output prefix"
  echo "  mq          minimum mapping quality requirements. Use 0 for all-to-dip, 1 for all-to-hap alignment."
  echo "  jid         OPTIONAL. Slurm jobid to wait."
  echo
  echo "  output      output files will be generated under dv_ONT_R9_MINMQ\${mq}_chr."
  exit -1
fi

ref=$1
out=$2
mq=$3
jid=$4
map=map-ont

if [[ $jid == "" ]]; then
  $tools/T2T-Polish/winnowmap/_submit.sh $ref $out $map
  jid=`cat filt.jid | tail -n1`
else
  echo "Skip submitting winnowmap."
fi

set -x
$tools/T2T-Polish/deepvariant/_submit_ont_r9_pepper_margin_dv.sh $out.pri.bam $ref $mq $jid

