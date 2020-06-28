#!/usr/bin/env bash

#SBATCH --array=0-31
#SBATCH --mem=10g
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --gres=lscratch:5
#SBATCH --partition=norm

######################################################################
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################

if [[ "$#" -lt 4 ]]; then
	echo "Usage: single_copy_filter.sh alignment target asm marker.meryl"
	echo
	echo "Filters alignment based on single-copy kmers"
	echo -e "\talignment: input bam or cram file, containing *only* target region (chr)."
	echo -e "\ttarget: target region (chr) to process"
	echo -e "\tasm: target fasta file"
	echo -e "\tmarker.meryl: meryl db containing marker kmers. This scripts generates chr specific markers."
	exit 0
fi


## Dependency: samtools, samToAlignment, meryl v1.0, subsetSamByKmers.py, samToErrorRate, SubFile, IGVTools

alignment=$1    # input bam / cram file 
target=$2	# chrX
asm=$3		# asm.fasta
single=$4	# single-copy.meryl
PREFIX=${alignment%%.*}	# Remove .*$ from $alignment

cores=$SLURM_CPUS_PER_TASK
if [ x$cores == "x" ]; then
	echo "Use 16 cores by default"
	cores=16
fi

#module load samtools

set -e

if [ ! -e $asm.fai ]; then
	echo "init.sh not executed"
	exit 0
fi
if [ ! -e $target.fa ]; then
	echo "init.sh not executed"
fi

if [ -e $target.markersandlength.cram ]; then
	echo "Already done"
	exit 0
else
	if [ ! -s $target.aligned.fasta ]; then

		# recompute variable EXPECTED before starting the for loop
		RID_TOTAL=`wc -l rid.list | awk '{print $1}'`
		echo "$RID_TOTAL"
		if [[ $RID_TOTAL -lt 100000 ]]; then
			EXPECTED=1
		else
			EXPECTED=$(($RID_TOTAL/100000))
			EXPECTED=$((EXPECTED+1))
		fi
		echo "Assuming $target.srt_id.bam was split into $EXPECTED files, with prefix split.$target.srt_id"

		# Sort by position -- this will be later distributed to job arrays
		for i in $(seq 1 $EXPECTED)
		do
			modulo_=$(($i%SLURM_ARRAY_TASK_COUNT))  #slurm job ids must be 0-based, e.g., --array=0-31
			if [[ $modulo_ -eq $SLURM_ARRAY_JOB_ID ]]; then
				i=$SLURM_ARRAY_TASK_ID
				split=split.$target.srt_id.$i
				echo "Processing file $split"
				cat $PREFIX.header > $split.sam
				cat $split | awk -v rid="a" '{if (rid!=$1) {print $1; rid=$1}}' > $split.tmp.list
				for rid in $(cat $split.tmp.list)
				do
					echo "Sort alignments from read $rid by position"
					cat $split | grep $rid | sort -nk2 >> $split.sam
				done
				echo -e "\nsamToAlignment $split.sam $asm"
				src/samToAlignment $split.sam $asm 2> $split.err \
					| awk '{if ($9 == 0) { print ">"$1"_"$5"_"$9"_"$10"_"$(NF-1); print $NF } else { print ">"$1"_"$5"_"$9"_"$12-$11"_"$(NF-1); print $NF}}' \
					> $split.aligned.fasta # && rm $split.sam $split.err $split.tmp.list
			fi
		done
		echo
	fi
fi
