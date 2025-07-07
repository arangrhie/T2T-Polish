#!/bin/bash

ref=$1
cpus=$SLURM_CPUS_PER_TASK

if [ -z $ref ]; then
	echo "Usage: ./bwa_index.sh <ref.fasta>"
	echo "No <ref.fasta> given. Exit."
	exit -1
fi

module load bwa
module load samtools

echo "
bwa index $ref"
bwa index $ref

if [[ ! -s $ref.fai ]]; then
  echo "samtools faidx $ref"
  samtools faidx -@$cpus $ref
fi

