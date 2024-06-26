#!/bin/bash

if ! [[ $# -eq 2 ]]; then
  echo "Usage: ./liftOverBed.sh in.bed chain"
  exit -1
fi

in=$1
in=${in/.bed/}
chain=$2

module load ucsc
set -x
liftOver ${in}.bed $chain ${in}.lifted.bed ${in}.lifted.unmapped.bed

