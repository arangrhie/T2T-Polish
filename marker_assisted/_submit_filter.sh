#!/bin/sh

SCRIPT=`cat SCRIPT`

if [[ "$#" -lt 5 ]]; then
	echo "Usage: ./submit_filter.sh in.bam target asm.fasta single.meryl array_len [len_filt]"
	echo -e "This script will submit the filtering and merging scripts"
	echo -e "\tin.bam: alignment file. .bam or .cram"
	echo -e "\ttarget: Chromosome or sequence id to extract"
	echo -e "\tasm.fasta: target.fasta or multi-fasta asm.fasta that contains the target sequence"
	echo -e "\tsingle.meryl: single-copy marker meryl db"
	echo -e "\tarray_len: to how many array jobs to submit"
	echo -e "\t[len_filt]: length filter. use LEN_FILT by default."
	echo 0
fi

bam=$1
target=$2
fa=$3
meryldb=$4
ARRAY_LEN=$5
len_filt=$6

mkdir -p logs

# computing requirements
cpus=4
mem=12g
partition=norm
walltime=1-0
path=`pwd`

# variables
name=$target.convert
args="$target $fa"
log=logs/$name.%A_%a.log
script=$SCRIPT/convert.sh
extra="--array=1-$ARRAY_LEN"

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > convert.jid

## Submit merge.sh
# wait until filt.sh finishes
cpus=24
mem=32g
jid=`cat convert.jid`
extra="--dependency=afterok:$jid"

name=$target.merge
args="$bam $target $fa $meryldb $len_filt"
log=logs/$name.%A_%a.log
script=$SCRIPT/merge.sh

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > merge.jid

