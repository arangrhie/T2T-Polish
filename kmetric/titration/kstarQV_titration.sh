#!/bin/bash

wait_file() {
  local file="$1"; shift

  until [ -f $file ] ; do sleep 300; done
  
}

seed=$RANDOM

mkdir -p exp${3}_rep${4}

seqtk sample -s$seed ${1} ${3} | gzip > exp${3}_rep${4}/fw.fastq.gz
seqtk sample -s$seed ${2} ${3} | gzip > exp${3}_rep${4}/rv.fastq.gz

cd exp${3}_rep${4}

ls *.fastq.gz > input.fofn

$VGP_PIPELINE/meryl2/_submit_meryl2_build.sh 21 input.fofn exp${3}_rep${4}

peak=$(meryl histogram exp${3}_rep${4}.meryl 2>/dev/null | awk '{if($2 > prev){dir="down"}else{dir="up"};if(prev_dir!=dir){counter+=1} prev_dir=dir;prev=$2;if (counter==4){print $1-1; exit}}')

wait_file exp${3}_rep${4}.meryl.hist

ln -s ../asm.meryl

sh ../qv.sh asm.meryl exp${3}_rep${4}.meryl $peak 32 > exp${3}_rep${4}.qv