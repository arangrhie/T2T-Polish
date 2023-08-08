#! /bin/bash

if [[ "$#" -lt 4 ]]; then
  echo "Usage: ./_submit_bwa_dv.sh ref.fasta fastq.map out sample [line_num]"
  echo "  ref.fasta  reference"
  echo "  fastq.map  space or tab teliminted for r1 and r2 in each line"
  echo "             example: /path/to/r1.fq.gz /path/to/r2.fq.gz"
  echo "  out        output prefix. used in bam"
  echo "  sample     sample name to appear in the final VCF"
  echo "  line_num   line_num to launch the mapping OPTIONAL"
  echo
  echo "  output VCF file will be generated as dv_mode.vcf.gz"
  echo
  echo "Not enough parameters given. Exit."
  exit -1
fi

ref=$1
fastq_map=$2
out=$3
sample=$4
mode="WGS" # no need to customize this, mode is always WGS

PIPELINE=$tools/T2T-Polish

sh $PIPELINE/bwa/_submit_bwa.sh $ref $fastq_map $out $5
wait_for=`tail -n1 mrg.jid`

sh $PIPELINE/deepvariant/_submit_deepvariant.sh $ref $out.dedup.bam $mode $sample $wait_for

