#! /bin/bash

if [[ "$#" -lt 5 ]] ; then
  echo "Usage: ./_submit_deepvariant_with_minqual.sh ref.fa alignment.bam mode sample minqual [jid_to_wait]"
  echo "  ref.fa          reference fasta file"
  echo "  alignment.bam   input alignment bam file"
  echo "  mode            WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  sample          name to be appear in the output VCF SAMPLE field"
  echo "  minqual         Recommending 0 for dip reference, -1 for haploid reference. 5 will be used for -1"
  echo "  output file will be generated as dv_[mode]_MQ[minqual].vcf.gz"
  exit -1
fi

REF=$(realpath $1)
BAM=$(realpath $2)
MODE=$3
SAMPLE=$4 # name to be appeared in the output VCF SAMPLE field
MINQUAL=$5 # optional
wait_for=$6 # optional
N_SHARD=12 # Lower to reduce racing condition

PIPELINE=$tools/T2T-Polish

if [[ $MINQUAL -eq -1 ]]; then
  MINQUAL=5
  echo "Using minqual=5 for haploid reference alignments."
fi

set -e
set -o pipefail

echo -e "$REF $BAM $MODE $SAMPLE"

# Check $MODE is valid
if [[ $MODE == "WGS" || $MODE == "PACBIO" || $MODE == "ONT_R104" || $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  echo "DeepVariant v1.6.1 in $MODE mode"
else
  echo "Unknown option $MODE"
fi

echo "Checking $REF.fai exists..."
if [[ ! -e $REF.fai ]]; then
  echo "No .fai found. Create one before launching with samtools faidx $REF"
  echo "Exit."
  exit -1
fi
echo

# Keep MODE and N_SHARD
echo $MODE > MODE_MQ$MINQUAL
echo $N_SHARD > N_SHARD_MQ$MINQUAL
echo $REF > REF_MQ$MINQUAL
echo $BAM > BAM_MQ$MINQUAL
echo $SAMPLE > SAMPLE_MQ$MINQUAL

mkdir -p logs

path=$PWD

extra=""

if [ ! -f deepvariant.step1.done ]; then
  echo "Step 1. make_examples"
  cpus=$N_SHARD
  mem=$(($cpus*4))g
  gres="lscratch:1000" # use /lscratch/ allow up to 1000 GB
  name=dv_step1_mq$MINQUAL
  script=$PIPELINE/deepvariant/step1_with_minqual.sh
  args="$SAMPLE $MINQUAL"
  partition=norm
  walltime=3-00:00:00
  log=logs/$name.%A.log

  if [[ ! -z $wait_for ]]; then
    extra="--dependency=afterok:$wait_for"
  fi

  set -x
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres  \
         --partition=$partition \
         -D $path \
         $extra --time=$walltime \
         --error=$log --output=$log $script $args > $name.jid
  set +x
  echo
  extra="--dependency=afterok:"`cat $name.jid`
else
  extra=""
fi
echo

if [ ! -f deepvariant.step2.done ]; then
  echo "Step 2. variant_call"

  cpus=12
  gres="gpu:1" # use /lscratch/ allow up to 50GB
  name=dv_step2_mq$MINQUAL
  script=$PIPELINE/deepvariant/step2_with_minqual.sh
  args="$MINQUAL"
  walltime=12:00:00
  partition=gpu
  log=logs/$name.%A.log
  
  echo
  set -x
  sbatch -J $name --cpus-per-task=$cpus --gres=$gres  \
         --partition=$partition \
         -D $path \
         $extra --time=$walltime \
         --error=$log --output=$log $script $args > $name.jid
  set +x
  extra="--dependency=afterok:"`cat $name.jid`
else
  extra=""
fi
echo

if [ ! -f deepvariant.step3.done ]; then
  echo "Step 3. postprocess_variants"
  cpus=12
  mem=120g
  name=dv_step3_mq$MINQUAL
  script=$PIPELINE/deepvariant/step3_with_minqual.sh
  args="$MINQUAL"
  partition=norm
  walltime=12:00:00
  log=logs/$name.%A.log

  echo
  set -x
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem  \
         --partition=$partition \
         -D $path \
         $extra --time=$walltime \
         --error=$log --output=$log $script $args
  set +x
fi
