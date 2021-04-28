#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./map.sh ref map opt [i]"
  echo "Map reads with winnowmap."
  echo "Requires input.fofn and repetitive_k15.txt"
  echo "  ref    : ref.fa"
  echo "  map    : Provide one from map-ont map-pb map-pb-clr"
  echo "  [i]    : fq file to proceed in line i of input.fofn"
  echo
  echo "  output : input fq prefix appended with .sam and an indexed .sort.bam"
  exit -1
fi

cpus=$(($SLURM_CPUS_PER_TASK-2))
if [[ -z $cpus ]]; then
  cpus=12
fi

i=$SLURM_ARRAY_TASK_ID
if [[ -z $i ]]; then
  i=$4
fi

ref=$1
map=$2
opt=$3

reads=`sed -n ${i}p input.fofn`

out=`basename $reads`
out=`echo $out | sed 's/.fasta.$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
out=`echo $out | sed 's/.fastq.$//g' | sed 's/.fastq$//g' | sed 's/.fq$//g' | sed 's/.fastq.gz$//g' | sed 's/.fq.gz$//g'`
out=$out.$i

if ! [[ -s $out.sam ]]; then
  echo "
  Align
  $tools/Winnowmap/bin/winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t$cpus $ref $reads > $out.sam
  "
  $tools/Winnowmap/bin/winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t$cpus $ref $reads > $out.sam
fi


if ! [[ -s $out.sort.bam ]]; then
  module load samtools

  echo "
  Sort
  samtools sort -@$cpus -m2G -T $out.tmp -O bam -o $out.sort.bam $out.sam"
  samtools sort -@$cpus -m2G -T $out.tmp -O bam -o $out.sort.bam $out.sam

  echo "
  Index"
  samtools index $out.sort.bam
fi

rm $out.sam

