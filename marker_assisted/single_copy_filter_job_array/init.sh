#!/usr/bin/env bash

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
	samtools faidx $asm
fi
if [ ! -e $target.fa ]; then
	samtools faidx $asm $target > $target.fa
	samtools faidx $target.fa
fi

if [ -e $target.markersandlength.cram ]; then
	echo "Already done"
	exit 0
else
	if [ ! -e $target.srt_id.bam ]; then
		echo "Get $target.srt_id.bam, sorted by read name"
		samtools sort -T $target.tmp -O bam -@$cores -m2G $alignment > $target.srt_id.bam

		# Get header
		samtools view -H $alignment > $PREFIX.header
		echo
	fi

	if [ ! -s $target.aligned.fasta ]; then
		echo "Split $target.srt_id.bam per 100000 reads"
		echo "# Get unique read ids"
		samtools view -O sam -@$cores $target.srt_id.bam | awk -v rid="a" '{if (rid!=$1) {print $1; rid=$1}}' > rid.list

		RID_TOTAL=`wc -l rid.list | awk '{print $1}'`
		echo "$RID_TOTAL"
		if [[ $RID_TOTAL -lt 100000 ]]; then
			echo "Total num. read ids are $RID_TOTAL : No need to split"
			samtools view -O sam -@$cores $target.srt_id.bam > split.$target.srt_id.1
			EXPECTED=1
		else
			EXPECTED=$(($RID_TOTAL/100000))
			for i in $(seq 1 $EXPECTED)
			do
				sed -n ${i}p rid.list >> rid.10kth	# 100000th read id
			done
			samtools view -@$cores $target.srt_id.bam | java -jar -Xmx4g src/txtSplitByFile.jar - rid.10kth 1 split.$target.srt_id
			EXPECTED=$((EXPECTED+1))
		fi
		echo "$target.srt_id.bam is split into $EXPECTED files, with prefix split.$target.srt_id"
		echo

		#For next steps, see main.sh and final.sh
	fi
fi
