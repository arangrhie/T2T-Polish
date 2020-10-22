#!/bin/bash

echo "Usage: ~/codes/sam2paf.sh in.bam out.paf"

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
  samtools view -h -@$cpu $1 | $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf - | awk \'{print \$1 .. \$16}\' - > $2"
  samtools view -h -@$cpu $1 | $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf - | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' - > $2
else
  echo "
  $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf $1 | awk \'{print \$1 .. \$16}\' - > $2"
  $tools/k8/k8 /usr/local/apps/minimap2/2.17/misc/paftools.js sam2paf $1 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' - > $2
fi
