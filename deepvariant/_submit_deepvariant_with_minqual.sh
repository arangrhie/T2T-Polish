#! /bin/bash

if [[ "$#" -lt 5 ]] ; then
  echo "Usage: ./_submit_deepvariant_with_minqual.sh ref.fa alignment.bam mode sample minqual [jid_to_wait]"
  echo "  mode   WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  output file will be generated as dv_mode.vcf.gz"
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
  mem=$(($cpus*3))g
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
