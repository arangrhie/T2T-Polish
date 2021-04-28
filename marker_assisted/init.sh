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
	echo "Usage: init.sh alignment target asm marker.meryl"
	echo
	echo "Filters alignment based on single-copy kmers"
	echo -e "\talignment: input bam or cram file, containing *only* target region (chr)."
	echo -e "\ttarget: target region (chr) to process"
	echo -e "\tasm: target fasta file"
	echo -e "\tmarker.meryl: meryl db containing marker kmers. This scripts generates chr specific markers."
	echo -e "\t[len_filt]: alignment length filter in kb. ex. 1 for 1kb"
	exit 0
fi


## Dependency: samtools, samToAlignment, meryl v1.0, subsetSamByKmers.py, samToErrorRate, SubFile, IGVTools

alignment=$1    # input bam / cram file 
target=$2	# chrX
asm=$3		# asm.fasta
meryldb=$4	# single-copy.meryl
len_filt=$5

cores=$SLURM_CPUS_PER_TASK
if [ x$cores == "x" ]; then
	echo "Use 16 cores by default"
	cores=16
fi

SCRIPT=`cat SCRIPT`
NUM_READS_PER_FILE=10000	# Increase or decrease depending on the read set
#NUM_READS_PER_FILE=100000	# Increase or decrease depending on the read set - for illumina reads

echo "Extract $target aslignments from $alignment, split by $NUM_READS_PER_FILE unique read ids"
echo

module load samtools

set -e

if [ ! -e $asm.fai ]; then
	samtools faidx $asm
fi
if [ ! -e $target.fa ]; then
	samtools faidx $asm $target > $target.fa
	samtools faidx $target.fa
fi

if [[ -z $len_filt ]]; then
	len_filt=`cat LEN_FILT`
fi

if [[ -z $len_filt ]]; then
	echo "LEN_FILT not found. Submit when running init.sh"
fi

if [ -e $target.markersandlength.cram ]; then
	echo "Already done"
	exit 0
else
	if [[ ! -s $target.bam ]]; then
		echo "Get $target.bam"
		samtools view -hb -@$cores $alignment $target > $target.bam
	fi
	if [[ ! -s $target.srt_id.bam ]]; then
		echo "Get $target.srt_id.bam, sorted by read name"
		samtools sort -n -T $target.tmp -O bam -@$cores -m1G $target.bam > $target.srt_id.bam

		# Get header
		samtools view -H $target.bam > $target.header
		echo
	fi

	if [ ! -s $target.rid.list ]; then
		echo "Get unique read ids"
		samtools view -O sam -@$cores $target.srt_id.bam | awk -v rid="a" '{if (rid!=$1) {print $1; rid=$1}}' > $target.rid.list
	fi

	if [ ! -s $target.RID_TOTAL ]; then
		RID_TOTAL=`wc -l $target.rid.list | awk '{print $1}'`
		echo $RID_TOTAL > $target.RID_TOTAL
	else
		RID_TOTAL=`cat $target.RID_TOTAL`
	fi
	
	echo "Total num. of unique reads aligned to $target: $target.$RID_TOTAL"

	EXPECTED=$(($RID_TOTAL/$NUM_READS_PER_FILE+1))
	FOUND=`ls split.$target.srt_id.* 2> /dev/null | wc -l`
	
	if [[ $FOUND -eq $EXPECTED ]]; then
		echo "*** Found $FOUND split.$target.srt_id.* files, matching expected $EXPECTED ***"
	else
		if [[ $EXPECTED -eq 1 ]]; then
			echo "$RID_TOTAL is less than $NUM_READS_PER_FILE : No need to split"
			samtools view -O sam -@$cores $target.srt_id.bam > split.$target.srt_id.1
		else
			echo "Collect every $NUM_READS_PER_FILE read id to $target.rid.ith"

			if [[ -s $target.rid.ith ]]; then
				echo "*** Found rid.ith. Removing it. ***"
				rm $target.rid.ith
			fi
			awk -v NUM_READS_PER_FILE=$NUM_READS_PER_FILE 'NR%NUM_READS_PER_FILE == 0 {print $0}' $target.rid.list >> $target.rid.ith

			echo "Split $target.srt_id.bam to $EXPECTED files"
			samtools view -@$cores $target.srt_id.bam | java -jar -Xmx4g $SCRIPT/src/txtSplitByFile.jar - $target.rid.ith 1 split.$target.srt_id
		fi
	fi
	echo

	echo "Submit _submit_filter.sh"
	echo "
	$SCRIPT/_submit_filter.sh $target.bam $target $asm $meryldb $EXPECTED $len_filt"
	$SCRIPT/_submit_filter.sh $target.bam $target $asm $meryldb $EXPECTED $len_filt
	echo
fi

echo "Clean up"
rm $target.srt_id.bam
