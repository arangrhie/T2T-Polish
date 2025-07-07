#! /bin/sh

if [[ "$#" -lt 6 ]]; then
  echo "Usage: ./_submit_mrg_dv.sh ref.fasta hifi.bam ilmn.bam out sample mq"
  echo "  ref.fasta  reference"
  echo "  hifi.bam   HiFi primary read alignment"
  echo "  ilmn.bam   Illumina read alignment. Use markdup.bam"
  echo "  out        output prefix for the merged bam"
  echo "  sample     sample name to appear in the final VCF"
  echo "  mq         minimum mapping quality requirements. Use positive values for MQ filters. Use -1 to follow pre-set options in DV."
  echo "             e.g. 0 for all-to-dip, -1 for all-to-hap alignments."
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
mq=$6

mode="HYBRID_PACBIO_ILLUMINA"

path=$PWD
cpus=12
mem=24g
partition=norm
walltime=12:00:00
script=$tools/T2T-Polish/deepvariant/merge_hybrid.sh
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
fi

name=sam2paf.$out
log=logs/$name.%A.log
script=$tools/T2T-Polish/coverage/sam2paf.sh
args="$out.bam $out.pri.paf"

wait_for=`tail -n1 mrg.jid`
extra="--dependency=afterok:$wait_for"

set -x
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \
       --partition=$partition \
       -D $path $extra --time=$walltime \
       --error=$log --output=$log $script $args
set +x

sh $tools/T2T-Polish/deepvariant/_submit_deepvariant_with_minqual.sh $ref $out.bam $mode $sample $mq $wait_for

