#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ast.sh <platform> [mean_cov]"
	exit -1
fi

platform=$1
mean_cov=$2
postfix=""

if [[ ! -z $mean_cov ]]; then
    mean_cov=`echo $mean_cov | awk '{printf "%.0f\n", $1*2.5}'`
    mean_cov="-M $mean_cov"
    postfix="_M"
fi

# Download and install from https://github.com/dfguan/asset/ . No third party dependency required other than zlib to run this coverage code
export asset=/data/Phillippy/tools/asset
pafs=`ls *.paf`

bed=${platform}${postfix}.bed

if [[ ! -s $bed ]]; then
echo "
	$asset/bin/ast_pb $mean_cov $pafs > $bed"
	$asset/bin/ast_pb ${mean_cov} $pafs > $bed && echo "Asset Finished!" || exit -1
fi

#echo "clean up"
#rm *.paf

module load bedtools

bedtools subtract -a ../asm.bed -b $bed | bedtools merge -d 100 -i - > $platform.low_high.bed

# Remove ends: remove low coverage region around ~1kb
bedtools subtract -a $platform.low_high.bed -b ../asm.ends.bed -A > $platform.low_high.trim1k.bed

# Get support: for getting reliable blocks. Will trim at the end when getting acc.reliable.bed
bedtools subtract -a ../asm.bed -b $platform.low_high.bed > $platform.support.bed

echo "# Rename the last pb.cov.wig accordingly and edit the track name"

