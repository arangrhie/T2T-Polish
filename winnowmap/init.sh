#!/bin/bash

if [[ -z $1 ]]; then
  echo "Usage: ./init.sh ref.fa"
  echo -e "\tCollect repetitive k=15 mers for winnowmap"
  exit -1
fi

ref=$1

meryl count k=15 $ref output merylDB  # add "compress" for homopolymer compression
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

