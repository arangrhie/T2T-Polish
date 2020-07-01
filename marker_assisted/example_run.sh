#!/bin/sh

# ./single_copy_filter.sh clr_test.bam chr20 t2t-chm13.20200602.fasta 20200602.single.meryl

# Submit to cluster
./_submit.sh clr.bam chr20 t2t-chm13.20200602.fasta 20200602.single.meryl 1
