#! /bin/bash

if [[ "$#" -lt 4 ]] ; then
  echo "Usage: ./_submit_deeptrio.sh ref.fa child.bam parent1.bam parent2.bam mode childsample parent1sample parent2sample [jid_to_wait]"
  echo "  mode   WGS | PACBIO | ONT_R104 | HYBRID_PACBIO_ILLUMINA"
  echo "  output file will be generated as dv_mode.vcf.gz"
  exit -1
fi

REF=$1
CHILDBAM=$2
PARENT1BAM=$3
PARENT2BAM=$4
MODE=$5
CHILDSAMPLE=$6 # name to be appeared in the output VCF SAMPLE field
PARENT1SAMPLE=$7 # name to be appeared in the first parent's output VCF SAMPLE field
PARENT2SAMPLE=$8 # name to be appeared in the second parent's output VCF SAMPLE field
wait_for=$9 # optional
N_SHARD=12

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
echo $CHILDBAM > CHILDBAM
echo $PARENT1BAM > PARENT1BAM
echo $PARENT2BAM > PARENT2BAM
echo $CHILDSAMPLE > CHILDSAMPLE
echo $PARENT1SAMPLE > PARENT1SAMPLE
echo $PARENT2SAMPLE > PARENT2SAMPLE

mkdir -p logs

for CHROM in `awk -F"\t" '{print $1}' $REF.fai`; do
    if [[ ! -d $CHROM\_$MODE\_output/examples ]]; then
        echo $CHROM " Step 1. make_examples"
        cpus=$N_SHARD
        mem=$(($cpus*2))g
        gres="lscratch:600" # use /lscratch/ allow up to 600 GB
        name=dv_step1
        script=$PIPELINE/deepvariant/deeptrio_step1.sh
        args="$REF $CHILDBAM $PARENT1BAM $PARENT2BAM $CHILDSAMPLE $PARENT1SAMPLE $PARENT2SAMPLE $CHROM"
        partition=norm
        walltime=12:00:00
        path=`pwd`
        log=logs/$name.%A.log
      
        if [[ ! -z $wait_for ]]; then
            extra="--dependency=afterok:$wait_for"
	else
	    extra=""
        fi
      
        echo
        set -x  
        sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres  \
               --partition=$partition \
               -D $path \
               $extra --time=$walltime \
               --error=$log --output=$log $script $args >> step1.jid
       set +x
        extra="--dependency=afterok:"`tail -1 step1.jid`
        
    else
        extra=""
    fi

    echo $CHROM " Step 2. variant_call"

    cpus=2
    mem=$(($cpus*2))g
    path=`pwd`
    gres="lscratch:50,gpu:p100:1" # use /lscratch/ allow up to 50GB
    name=dv_step2
    script=$PIPELINE/deepvariant/deeptrio_step2.sh
    walltime=4:00:00
    args="$CHROM"
    partition=gpu
    log=logs/$name.%A.log
    
    echo "
    sbatch -J $name --cpus-per-task=$cpus --mem=$mem --gres=$gres \\
           --partition=$partition \\
           -D $path \\
           --time=$walltime $extra \\
           --error=$log --output=$log $script $args >> step2.jid
    "
    
    sbatch -J $name --cpus-per-task=$cpus --gres=$gres --mem=$mem \
           --partition=$partition \
           -D $path \
           --time=$walltime $extra \
           --error=$log --output=$log $script $args >> step2.jid
    
    echo "Step 3. postprocess_variants"
    cpus=8
    mem=48g
    name=dv_step3
    script=$PIPELINE/deepvariant/deeptrio_step3.sh
    args="$CHROM"
    partition=norm
    walltime=12:00:00
    extra="--dependency=afterok:"`tail -1 step2.jid`
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
           --error=$log --output=$log $script $args >> step3.jid
done
