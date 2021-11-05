#!/bin/bash

printf "merfin -sequence ${1} -seqmers ${2} -readmers ${3} -peak ${4} -vcf ${5} -output merfin/${6} -memory 200 -vmer\n\n"
merfin -sequence ${1} -seqmers ${2} -readmers ${3} -peak ${4} -vcf ${5} -output ${6} -memory1 30 -memory2 200 -vmer
