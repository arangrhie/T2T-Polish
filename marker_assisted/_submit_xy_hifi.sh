#!/bin/bash

OUT=hg002XYv2.7_hifi.pri
BAM="$OUT.bam"
FA_X="$projects/HG002_T2T/assembly/v2.7.X.fasta"
FA_Y="$projects/HG002_T2T/assembly/v2.7.Y.fasta"
MAT_MER="$projects/HG002_T2T/markers/v2.7.k21.marker.mat.meryl"
PAT_MER="$projects/HG002_T2T/markers/v2.7.k21.marker.pat.meryl"
LEN=10

ln -s ../$BAM
ln -s ../$BAM.bai
ln -s $MAT_MER
ln -s $PAT_MER

set -x
sh ~/codes/_submit_norm.sh 24 32g chrX $tools/T2T-Polish/marker_assisted/filter_by_marker_nosplit.sh "$BAM chrX $FA_X $MAT_MER $LEN" | tail -n1 | awk '{print $NF}' >  filter.jids
sh ~/codes/_submit_norm.sh 24 32g chrY $tools/T2T-Polish/marker_assisted/filter_by_marker_nosplit.sh "$BAM chrY $FA_Y $PAT_MER $LEN" | tail -n1 | awk '{print $NF}' >> filter.jids

jids=`cat filter.jids | tr '\n' ',' | sed 's/,$//g'`
sh ~/codes/_submit_norm.sh 12 24g aggregate $tools/T2T-Polish/marker_assisted/aggregate.sh $OUT "--dependency=afterok:$jids"
