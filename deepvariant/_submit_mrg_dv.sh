#! /bin/bash

if [[ "$#" -lt 5 ]]; then
  echo "Usage: ./_submit_mrg_dv.sh ref.fasta hifi.bam ilmn.bam out sample [line_num]"
  echo "  ref.fasta  reference"
  echo "  hifi.bam   HiFi primary read alignment"
  echo "  ilmn.bam   Illumina read alignment. Use markdup.bam"
  echo "  out        output prefix for the merged bam"
  echo "  sample     sample name to appear in the final VCF"
  echo "  line_num   line_num to launch the mapping OPTIONAL"
  echo
  echo "  output VCF file will be generated as dv_mode.vcf.gz"
  echo
  echo "Not enough parameters given. Exit."
  exit -1
fi

ref=$1
hifi_bam=$2
ilmn_bam=$3
out=$4
sample=$5
mode="HYBRID_PACBIO_ILLUMINA"

PIPELINE=$tools/T2T-Polish

path=$PWD
cpus=12
mem=24g
partition=norm
walltime=12:00:00
script=$PIPELINE/deepvariant/merge_hybrid.sh
args="$hifi_bam $ilmn_bam $out.bam"
name=mrg.$out
log=logs/$name.%A.log

set -o pipefail
set -e

if ! [[ -s $out.bam ]]; then

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem  \
       --partition=$partition \
       -D $path \
       $extra --time=$walltime \
       --error=$log --output=$log $script $args > mrg.jid
set +x

wait_for=`tail -n1 mrg.jid`
fi

name=sam2paf.$out
log=logs/$name.%A.log
script=$PIPELINE/coverage/sam2paf.sh
args="$out.bam $out.pri.paf"

jid=`tail -n1 filt.jid`
extra="--dependency=afterok:$wait_for"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \
       --partition=$partition \
       -D $path $extra --time=$walltime \
       --error=$log --output=$log $script $args
set +x

sh $PIPELINE/deepvariant/_submit_deepvariant.sh $ref $out.bam $mode $sample $wait_for

