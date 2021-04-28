#!/bin/sh

# Submit to cluster
for i in $(seq 1 22) X
do
  ./_submit.sh clr.bam chr$i t2t-chm13.20200602.fasta 20200602.single.meryl 1
done

# or run init.sh directly
./init.sh clr.bam chr20 t2t-chm13.20200602.fasta 20200602.single.meryl 1
