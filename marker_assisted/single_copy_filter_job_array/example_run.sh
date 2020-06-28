#!/bin/sh

./init.sh clr_test.bam chr20 t2t-chm13.20200602.fasta 20200602.single.meryl
sbatch main.sh clr_test.bam chr20 t2t-chm13.20200602.fasta 20200602.single.meryl
./final.sh clr_test.bam chr20 t2t-chm13.20200602.fasta 20200602.single.meryl
