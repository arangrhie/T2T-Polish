#!/bin/bash

if [[ "$#" -lt  2 ]]; then
  echo "Usage: ./merge.sh out-prefix in.list"
  echo "Merge bams from in.list to out-prefix.bam"
  exit -1
fi

cpu=$SLURM_CPUS_PER_TASK
out=$1
lst=$2

if [ -z $out ]; then
	exit -1
fi

module load samtools

set -e
set -o pipefail

bams=`cat $lst   | tr '\n' ' '`
csis=`echo $bams | sed 's/\.bam$/.bam.csi/g'`

num_bams=`wc -l $lst | awk '{print $1}'`

if [[ "$num_bams" -eq 1 ]]; then
	echo "Only 1 bam provided. Skipping Merging."
  mv $bams $out.dedup.pri.bam
  mv $csis $out.dedup.pri.bam.csi
  touch bwa.done
	exit 0
fi

echo "Merge $bams"

set -x
samtools merge -@ $cpu -O bam -b $lst $out.dedup.pri.bam
samtools index -@ $cpu $out.dedup.pri.bam
set +x

echo "
Clean up"
rm $bams $csis
touch bwa.done