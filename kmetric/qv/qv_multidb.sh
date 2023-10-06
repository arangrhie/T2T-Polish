#!/bin/sh

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./qv.sh asm.meryl read1.meryl peak1 read2.meryl peak2 threads"
	echo
	echo "  asm.meryl     meryl db of the assembly"
	echo "  read1.meryl   meryl db of the first read set"
        echo "  read2.meryl   meryl db of the second read set"
	echo "  peak1         haploid peak to normalize kmer multiplicity of the read.meryl"
        echo "  peak2         haploid peak to normalize kmer multiplicity of the read.meryl"
	echo "  threads       number of threads to use"
	echo
	echo "Output (stdout) format: asm.meryl <tab> excessive <tab> total <tab> qv <tab> error"
	echo "    asm.meryl   name of the provided meryl db of the assembly"
	echo "    excessive   excessive number of kmers found in the assembly"
	echo "    total       total number of kmers in the assembly"
	echo "    qv          estimated base level quality value"
	echo "    error       error rate used to calculate the qv"
	echo
	echo "Aug. 12, 2020"
	echo "Giulio Formenti gformenti@rockefeller.edu"
	echo "Arang  Rhie     arrhie@gmail.com"
	echo
	exit -1;
fi

asm=$1
red1=$2
p1=$3
red2=$4
p2=$5
threads=$6
if [ -z $threads ]; then
  threads=8
fi

meryl threads=$threads intersect $red1 $asm output ${red1}_and_$asm
meryl threads=$threads intersect $red2 $asm output ${red2}_and_$asm
meryl threads=$threads divide-round $p1 ${red1}_and_$asm output kr1.meryl
meryl threads=$threads divide-round $p2 ${red2}_and_$asm output kr2.meryl

meryl intersect-max kr1.meryl kr2.meryl output kr.meryl

meryl threads=$threads subtract $asm kr.meryl output ka-extra.meryl

k=`meryl threads=$threads print $asm | head -n 2 | tail -n 1 | awk '{print length($1)}'`

# echo "# QV statistics for $asm"
ASM_ONLY=`meryl statistics ka-extra.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
TOTAL=`meryl statistics $asm  | head -n4 | tail -n1 | awk '{print $2}'`
ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR"
