#!/bin/sh

SCRIPT=$tools/T2T-Polish/marker_assisted
echo $SCRIPT > SCRIPT

if [[ "$#" -lt 5 ]]; then
	echo "Usage: ./submit.sh in.bam target asm.fasta marker.meryl len_filt"
	echo -e "\tin.bam: alignment file. .bam or .cram"
	echo -e "\ttarget: chromosome or sequence id to extract"
	echo -e "\tasm.fasta: target.fasta or multi-fasta asm.fasta which contains the target sequence"
	echo -e "\tmarker.meryl: marker meryl db (ex. single-copy kmer db)"
	echo -e "\tlen_filt: length filter, in kb (INTEGER). ex. 1 for 1kb"
	exit 0
fi

bam=$1
target=$2
fa=$3
meryldb=$4
len_filt=$5
echo $len_filt > LEN_FILT

cpus=24
mem=48g
name=$target.init
script=$SCRIPT/init.sh

args="$bam $target $fa $meryldb"
partition=norm
walltime=1-0
path=`pwd`

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args

