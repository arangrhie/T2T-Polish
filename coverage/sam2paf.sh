#!/bin/bash

echo "Usage: ./sam2paf.sh in.bam out.paf [name]"

if [[ "$#" -lt 2 ]]; then
  echo "No input and output file provied."
  exit -1
fi

cpu=$SLURM_CPUS_PER_TASK

module load minimap2/2.17

# Download k8 under $tools/k8 (https://github.com/lh3/minimap2/tree/master/misc/README.md)
# cd $tools/k8
# curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
# ln -s k8-`uname -s` k8

if [[ $1 =~ .bam$ ]] || [[ $1 =~ .cram$ ]]; then
  echo "
  samtools view -h -@$cpu $1 | $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf - | cut -f1-16 - > $2"
  samtools view -h -@$cpu $1 | $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf - | cut -f1-16 - > $2
else
  ## Add -L for long cs tag form
  echo "
  $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf $1 | cut -f1-16 - > $2"
  $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf $1 | cut -f1-16 - > $2

fi

out_prefix=`echo $2 | sed 's/.paf//g'`
name=$3

if [[ -z $name ]]; then
  echo "No track name provided. Stopping here."
  exit 0
fi

echo "
java -jar -Xmx4g $tools/T2T-Polish/paf_util/pafToCovWig.jar $2 "$name" 1024 > $out_prefix.cov.wig"
java -jar -Xmx4g $tools/T2T-Polish/paf_util/pafToCovWig.jar $2 "$name" 1024 > $out_prefix.cov.wig
