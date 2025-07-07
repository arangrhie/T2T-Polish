#!/bin/bash

ref=$1
fastq_map=$2
out=$3
wd=$PWD

if [[ -z $ref || -z $fastq_map || -z $out ]]; then
	echo "Usage: ./bwa.sh <ref.fasta> <fastq_map> <out>"
	echo "No <ref.fasta> found. Exit."
	exit -1
fi

module load bwa
module load samtools # v1.21+

cpu=$SLURM_CPUS_PER_TASK
idx=$SLURM_ARRAY_TASK_ID
tmp=/lscratch/${SLURM_JOB_ID}

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
  # try using $4
  idx=$4
fi

if [[ -z "$idx" ]]; then
  echo "No idx found"
  exit -1
fi
out=$out.$idx

line=`sed -n ${idx}p $fastq_map`
r1=`echo $line | awk '{print $1}'`
r2=`echo $line | awk '{print $2}'`

set -o pipefail
set -e

# $out.bam already exists
if [[ -s $out.bam && ! -s $out.dedup.bam ]]; then
  echo "Found $out.bam, but no $out.dedup.bam"
  echo "Processing $out.bam without re-running bwa steps"

set -x
  samtools markdup -r -@$cpu $out.bam $out.dedup.bam
  samtools index $out.dedup.bam
set +x
  echo "Done!"
  exit 0
fi

set -x

bwa mem -t $cpu $ref $r1 $r2 > $tmp/$out.sam

samtools fixmate -m -@$cpu $tmp/$out.sam $tmp/$out.fix.bam && rm $tmp/$out.sam

samtools sort -@$cpu -O bam -o $tmp/$out.bam -T $tmp/$out.tmp $tmp/$out.fix.bam && rm $tmp/$out.fix.bam
samtools index -@$cpu $tmp/$out.bam

samtools markdup -r -@$cpu $tmp/$out.bam $tmp/$out.dedup.bam && rm $tmp/$out.bam
samtools index -@$cpu $tmp/$out.dedup.bam

# Filter out secondary alignments
samtools view -@$cpu -F0x100 -hb --write-index -o $tmp/$out.dedup.pri.bam $tmp/$out.dedup.bam &&
mv $tmp/$out.dedup.pri.* ${wd}/ || exit -1
rm -rf $tmp/*
touch $out.done

set +x

echo "Done!"
