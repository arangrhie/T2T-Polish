#! /bin/bash

if [[ "$#" -lt 4 ]] ; then
  echo "Usage: ./_submit_deepvariant.sh ref.fa alignment.bam mode sample [jid_to_wait]"
  echo "  mode   WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  output file will be generated as dv_mode.vcf.gz"
  exit -1
fi

REF=$1
BAM=$2
MODE=$3
SAMPLE=$4 # name to be appeared in the output VCF SAMPLE field
wait_for=$5 # optional
N_SHARD=48

PIPELINE=$tools/T2T-Polish

# Check $MODE is valid
if [[ $MODE == "WGS" || $MODE == "PACBIO" || $MODE == "ONT_R104" || $MODE == "HYBRID_PACBIO_ILLUMINA" ]]; then
  echo "DeepVariant v1.5 in $MODE mode"
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
echo $MODE > MODE
echo $N_SHARD > N_SHARD
echo $REF > REF
echo $BAM > BAM
echo $SAMPLE > SAMPLE

mkdir -p logs

if [[ ! -d dv_$MODE/examples ]]; then
  echo "Step 1. make_examples"
  cpus=$N_SHARD
  mem=$(($cpus*3))g
  gres="lscratch:1000" # use /lscratch/ allow up to 1000 GB
  name=dv_step1
  script=$PIPELINE/deepvariant/step1.sh
  args="$REF $BAM $SAMPLE"
  partition=norm
  walltime=12:00:00
  path=`pwd`
  log=logs/$name.%A.log

  if [[ ! -z $wait_for ]]; then
    extra="--dependency=afterok:$wait_for"
  fi

  echo "
  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres \\
         --partition=$partition \\
         -D $path \\
         $extra --time=$walltime \\
         --error=$log --output=$log $script $args
  "

  sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres  \
         --partition=$partition \
         -D $path \
         $extra --time=$walltime \
         --error=$log --output=$log $script $args > step1.jid
  extra="--dependency=afterok:"`cat step1.jid`
else
  extra=""
fi

echo "Step 2. variant_call"

cpus=12
gres="lscratch:50,gpu:p100:4" # use /lscratch/ allow up to 50GB
name=dv_step2
script=$PIPELINE/deepvariant/step2.sh
args=""
partition=gpu
log=logs/$name.%A.log

echo "
sbatch -J $name --cpus-per-task=$cpus --gres=$gres \\
       --partition=$partition \\
       -D $path \\
       $extra --time=$walltime \\
       --error=$log --output=$log $script $args > step2.jid
"

sbatch -J $name --cpus-per-task=$cpus --gres=$gres  \
       --partition=$partition \
       -D $path \
       $extra --time=$walltime \
       --error=$log --output=$log $script $args > step2.jid

echo "Step 3. postprocess_variants"
cpus=8
mem=48g
name=dv_step3
script=$PIPELINE/deepvariant/step3.sh
args=""
partition=norm
walltime=12:00:00
extra="--dependency=afterok:"`cat step2.jid`
log=logs/$name.%A.log

echo "
sbatch -J $name --cpus-per-task=$cpus --mem=$mem \\
       --partition=$partition \\
       -D $path \\
       $extra --time=$walltime \\
       --error=$log --output=$log $script $args
"
sbatch -J $name --cpus-per-task=$cpus --mem=$mem  \
       --partition=$partition \
       -D $path \
       $extra --time=$walltime \
       --error=$log --output=$log $script $args

