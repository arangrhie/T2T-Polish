#! /bin/bash

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./_submit_bwa_dv_illumina.sh ref out mq [line_num]"
  echo "  ref        reference"
  echo "  out        output prefix. used in bam"
  echo "  mq         minimum mapping quality requirements. Use positive values for MQ filters. Use -1 to follow pre-set options in DV."
  echo "             e.g. 0 for all-to-dip, -1 for all-to-hap alignments."
  echo "  line_num   line_num to launch the mapping OPTIONAL"
  echo
  echo "  output VCF file will be generated as dv_mode.vcf.gz"
  echo
  echo "Assumes input.fofn in the same path. REQUIRED."
  echo "  input.fofn space or tab teliminted for r1 and r2 in each line"
  echo "             example: /path/to/r1.fq.gz /path/to/r2.fq.gz"
  echo "Not enough parameters given. Exit."
  exit -1
fi

ref=$1
out=$2
mq=$3
line_num=$4
fastq_map=input.fofn
mode="WGS"      # no need to customize this, mode is always WGS
sample=illumina # SAMPLE name in the output VCF

set -x
sh $tools/T2T-Polish/bwa/_submit_bwa.sh $ref $fastq_map $out $line_num
wait_for=`tail -n1 mrg.jid`

if [[ "$mq" -eq -1 ]]; then
  sh $tools/T2T-Polish/deepvariant/_submit_deepvariant.sh $ref $out.dedup.bam $mode $sample $wait_for
else
  sh $tools/T2T-Polish/deepvariant/_submit_deepvariant_with_minqual.sh $ref $out.dedup.bam $mode $sample $mq $wait_for
fi

