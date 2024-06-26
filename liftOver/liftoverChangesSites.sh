#!/bin/bash

if [[ -z $2 ]]; then
	echo "Usage: ./liftoverChangesSites.sh in.vcf chain"
	echo "  output: in.bed, in.lifted.bed, in.lifted.unmapped.bed, in.lifted.srt.bed"
  echo "  in.bed : SUB at same POS + length(REF)"
  echo "           DEL at same POS"
  echo "           INS at same POS"
	exit -1
fi

in=$1
name=`echo $in | sed 's/.vcf.gz$//g' | sed 's/.vcf$//g' | sed 's/.bcf$//g'`

if [[ $in == *.gz ]]; then
  bcftools view -Ov $in > $name.vcf
fi

chain=$2

# DEL: vlen > 0
# SUB: vlen = 0
# INS: vlen < 0
# Start = $2-1 (POS-1) to get 0-base
# End   = $2 (POS) + length(REF) // SUB
#       = $2 (POS)               // DEL
#       = $2 (POS)               // INS

awk -F "\t" '$1 !~ /#/ {  \
  REF=length($4); \
  ALT=length($5); \
  vlen=REF-ALT; \
  if (vlen == 0) { note = "SUB"; start=$2-1; end=$2+REF; vlen=REF } \
    else if (vlen > 0) { note="DEL"; start=$2-1; end=$2; } \
    else { note = "INS"; start=$2-1; end=$2; vlen*=-1 } \
  {print $1"\t"start"\t"end"\t"note":"vlen":"REF":"ALT}}' \
  $name.vcf > $name.bed

module load ucsc
module load bedtools

set -e

# This step has some bugs around large indels. let's create a dummy and correct it
# since we know the ALT length
# liftOver ${name}.bed $chain ${name}.lifted.bed ${name}.lifted.unmapped.bed

# create a dummy for regions,
# this dummy has START=POS-2 or END=POS+REF+1 in bed format; all in length 1
# with ALT length in 8th column
cat $name.vcf | awk -F"\t" '$1 !~ /^#/ { \
  REF=length($4); \
  ALT=length($5); \
  POS=$2;         \
  vlen=REF-ALT;   \
  if (vlen == 0) { note = "SUB"; vlen=REF } \
    else if (vlen > 0) { note="DEL" } \
    else { note = "INS"; vlen*=-1 } \
  if (POS > 2) {print $1"\t"POS-2"\t"POS-1"\t"ALT":"note":"vlen":"REF":"ALT} \
  else         {print $1"\t"POS+REF"\t"POS+REF+1"\t0:"note":"vlen":"REF":"ALT} }' \
  > dummy.bed

liftOver dummy.bed $chain dummy.lifted.bed dummy.unlifted.bed

if [[ -s dummy.unlifted.bed ]]; then
  #  if dummy has problems
  echo "dummy has problems. try conventional liftOver for this case"

# shift it back to its original form 
  cat dummy.unlifted.bed | tr ':' '\t' |\
    awk -F"\t" '$1 !~ /^#/ { \
      if ($4>0) {print $1"\t"$3"\t"$3+$7"\t"$5":"$6":"$7":"$8} \
      else      {print $1"\t"$2-$7"\t"$2"\t"$5":"$6":"$7":"$8} \
    }' >> dummy2.bed

  liftOver dummy2.bed $chain dummy2.lifted.bed dummy2.unlifted.bed

  if [[ -s dummy2.unlifted.bed ]]; then
    echo "liftOver of the original is problematic again"
    cat dummy2.unlifted.bed
    exit -1
  fi

  cat dummy2.lifted.bed > $name.lifted.bed
fi

# add lifted regions, correct their START end END accordingly
cat dummy.lifted.bed | tr ':' '\t' |\
  awk '{ if ($4>0) {print $1"\t"$3"\t"$3+$8"\t"$5":"$6":"$7":"$8} \
         else      {print $1"\t0\t"$2"\t"$5":"$6":"$7":"$8}
  }'>> $name.lifted.bed

cat $name.lifted.bed | bedtools sort -i - > $name.lifted.srt.bed
mv $name.lifted.srt.bed $name.lifted.bed

rm dummy*

echo "Done!"
echo "Output written into $name.lifted.bed"

