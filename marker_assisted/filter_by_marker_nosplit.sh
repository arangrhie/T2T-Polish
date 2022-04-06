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

if [[ "$#" -lt 5 ]]; then
  echo "Usage: filter_by_marker_nosplit.sh alignment target asm marker.meryl len_filt"
  echo
  echo "Filters alignment based on single-copy kmers"
  echo -e "\talignment: input bam or cram file"
  echo -e "\ttarget: target region (chr) to process"
  echo -e "\tasm: entire fasta file"
  echo -e "\tmarker.meryl: meryl db containing marker kmers. This scripts generates chr specific markers."
  echo -e "\tlen_filt: alignment length filter in kb. ex. 1 for 1kb"
  exit 0
fi


## Dependency: samtools, samToAlignment, meryl v1.3, subsetSamByKmers.py, samToErrorRate, SubFile

set -e
set -x

alignment=$1  # input bam / cram file 
target=$2     # chrX
asm=$3        # asm.fasta
meryldb=$4    # single-copy.meryl
len_filt=$5

if [[ -z $len_filt ]]; then
  echo "No len_filt provided. Exit."
  exit -1
fi

cores=$SLURM_CPUS_PER_TASK
if [ x$cores == "x" ]; then
  echo "Use 16 cores by default"
  cores=16
fi

SCRIPT=$tools/T2T-Polish/marker_assisted
echo $SCRIPT > SCRIPT

echo $len_filt > LEN_FILT

if [ ! -e $asm.fai ]; then
  samtools faidx $asm
fi
if [ ! -e $target.fa ]; then
  samtools faidx $asm $target > $target.fa
  samtools faidx $target.fa
fi

if [ -e $target.markersandlength.bam ]; then
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

  if [[ ! -s $target.aligned.fasta ]]; then
    echo "Sort by flag, convert to extended fa"
    cat $target.header > $target.tmp.sam
    samtools view -@$cores $target.srt_id.bam | java -jar -Xmx12g $SCRIPT/src/samSortByFlag.jar - >> $target.tmp.sam
    $SCRIPT/src/samToAlignment $target.tmp.sam $target.fa 2> $target.tmp.err \
      | awk '{if ($9 == 0) { print ">"$1"_"$5"_"$9"_"$10"_"$(NF-1); print $NF } else { print ">"$1"_"$5"_"$9"_"$12-$11"_"$(NF-1); print $NF}}' \
      > $target.aligned.fasta

    # echo "Clean up $target.srt_id.bam $(ls $target.tmp.*)"
    rm $target.srt_id.bam $target.tmp.*
  fi
  
  if [ ! -s $target.alignment.posCount ]; then
    echo "Prepare $target.single.meryl"
    meryl count k=21 $target.fa output $target.meryl
    meryl intersect $target.meryl $meryldb output $target.single.meryl
    meryl-lookup -existence -sequence $target.aligned.fasta -mers $target.single.meryl | awk '$NF>0 {print $1"\t"$NF}' > $target.alignment.posCount
    echo
  fi

  echo "Subset by markers"
  if [ ! -s $target.markers.bam ]; then
    samtools view -h -O sam -@$cores $target.bam > $target.sam
    echo "python $SCRIPT/src/subsetSamByKmers.py $target.alignment.posCount $target.sam > $target.markers.sam"
    python $SCRIPT/src/subsetSamByKmers.py $target.alignment.posCount $target.sam > $target.markers.sam
    samtools view -@$cores -hb -o $target.markers.bam $target.markers.sam
  fi
  
  echo
  echo "# hard filter alignments < $len_filt kb to $target.filtered.sam"
  if [ -s $target.filtered.sam ]; then
    echo "*** Found $target.filtered.sam. Skipping this step. ***"
  else
    cat $target.header > $target.filtered.sam
    samtools view -@$cores $target.markers.sam |awk -v l=$len_filt '{if (length($10) >= l) print $0}' >> $target.filtered.sam
  fi

  echo
  echo "# also filter by alignment legnth and 75% idy"
  if [ -s $target.filteredList ]; then
    echo "*** Found $target.filteredList. Skipping this step. ***"
  else
    $SCRIPT/src/samToErrorRate $target.filtered.sam $asm \
      | python3 $SCRIPT/src/lenFiltUniq.py $len_filt 75 \
      > $target.filteredList
  fi

  cat $target.header > $target.tmp.sam
  samtools view -@$cores $target.filtered.sam > $target.tmp2.sam
  java -cp $SCRIPT/src/ SubFile $target.filteredList $target.tmp2.sam >> $target.tmp.sam
  mv $target.tmp.sam $target.filtered.sam
  rm $target.markers.sam $target.tmp2.sam
fi

echo "
# generate $target.markersandlength.bam"
samtools view -hb -@$cores -o $target.markersandlength.bam $target.filtered.sam

echo "
# Index"
samtools index $target.markersandlength.bam

echo "
# Cleanup"
rm -r $target.filtered.sam $target.sam $target.bam $target.align* $target.filteredList $target.meryl $target.single.meryl $target.fa $target.fa.fai $target.header

