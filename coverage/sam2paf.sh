#!/bin/bash

echo "Usage: ./sam2paf.sh in.bam out.paf [name]"

if [[ "$#" -lt 2 ]]; then
  echo "No input and output file provied."
  exit -1
fi

cpu=$SLURM_CPUS_PER_TASK

module load minimap2/2.26
module load samtools

# Download k8 under $tools/k8 (https://github.com/lh3/minimap2/tree/master/misc/README.md)
# cd $tools/k8
# curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
# ln -s k8-`uname -s` k8

mm2_path=`which minimap2`
mm2_path=`dirname $mm2_path`

set -x
if [[ $1 =~ .bam$ ]] || [[ $1 =~ .cram$ ]]; then
  samtools view -h -@$cpu $1 | $tools/k8/k8 $mm2_path/misc/paftools.js sam2paf - | cut -f1-16 - > $2
else
  ## Add -L for long cs tag form
  $tools/k8/k8 $mm2_path/misc/paftools.js sam2paf $1 | cut -f1-16 - > $2

fi
set +x

out_prefix=`echo $2 | sed 's/.paf//g'`
name=$3

if [[ -z $name ]]; then
  echo "No track name provided. Stopping here."
  exit 0
fi

set -x
java -jar -Xmx4g $tools/T2T-Polish/paf_util/pafToCovClippedWig.jar $2 "$name" 1024 $out_prefix
