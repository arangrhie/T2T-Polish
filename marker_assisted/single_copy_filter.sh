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

echo "*** DEPRECATED. USE _submit.sh instead. ***"
echo


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
      samtools sort -n -T $target.tmp -O bam -@$cores -m2G $alignment > $target.srt_id.bam

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

	      echo "Collect every $NUM_READS_PER_FILE read id to rid.ith"

	      if [[ -s rid.ith ]]; then
	         echo "*** Found rid.ith. Removing it. ***"
	         rm rid.ith
	      fi
	      awk -v NUM_READS_PER_FILE=$NUM_READS_PER_FILE 'NR%NUM_READS_PER_FILE == 0 {print $0}' rid.list >> rid.ith

	      samtools view -@$cores $target.srt_id.bam | java -jar -Xmx4g src/txtSplitByFile.jar - rid.ith 1 split.$target.srt_id
	      EXPECTED=$((EXPECTED+1))
	   fi
	   echo "$target.srt_id.bam is splitted into $EXPECTED files, with prefix split.$target.srt_id"
	   echo

	   # Sort by position -- this will be later distributed to job arrays
	   for i in $(seq 1 $EXPECTED)
	   do
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
	   done
	   echo

	   echo "Merge to $target.aligned.fasta"
	   cat split.$target.srt_id.*.aligned.fasta > $target.aligned.fasta
	   echo
   fi

   if [ ! -s $target.alignment.posCount ]; then
      echo "Prepare $target.single.meryl"
      meryl count k=21 $target.fa output $target.meryl
      meryl intersect $target.meryl $single output $target.single.meryl
      meryl-lookup -existence -memory 12 -sequence $target.aligned.fasta -mers $target.single.meryl | awk '$NF>0 {print $1"\t"$NF}' > $target.alignment.posCount
      echo
   fi

   echo "Subset by markers"
   samtools view -h -O sam -@$cores $alignment > $target.sam
   echo "python src/subsetSamByKmers.py $target.alignment.posCount $target.sam > $target.markers.sam"
   python src/subsetSamByKmers.py $target.alignment.posCount $target.sam > $target.markers.sam
   samtools sort -@$cores -O cram -o $target.markers.cram -T $target.$PREFIX.tmp --reference=$asm $target.markers.sam
   #samtools index $target.markers.cram
   $tools/IGVTools/igvtools count $target.markers.cram $target.markers.tdf $asm.fai
   echo

   echo "# filter alignments <50kb to $target.filtered.sam"
   cat $PREFIX.header > $target.filtered.sam
   cat $target.markers.sam |grep -v "^@" |awk '{if (length($10) >= 50000) print $0}' >> $target.filtered.sam

   # also filter by alignment legnth
   src/samToErrorRate $target.filtered.sam $asm \
      | awk '{if ($9 == 0) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t1\t"$12-$11"\t"$12-$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' \
      | awk '{if ($3 <= -25000 && $4 >= 75) print $target}'| awk -v rid="" '{if (rid!=$1) {print $1; rid=$1}}' > $target.$PREFIX.filteredList
   cat $PREFIX.header > $target.tmp.sam
   cat $target.filtered.sam |grep -v "^@" > $target.tmp2.sam
   java src/SubFile $target.$PREFIX.filteredList $target.tmp2.sam >> $target.tmp.sam
   mv $target.tmp.sam $target.filtered.sam
fi

samtools sort -@${cores} -O cram -o $target.markersandlength.cram -T $target.$PREFIX.tmp --reference=$asm $target.filtered.sam
#samtools index $target.markersandlength.cram
$tools/IGVTools/igvtools count $target.markersandlength.cram $target.markersandlength.tdf $asm.fai

