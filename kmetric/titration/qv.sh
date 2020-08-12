#!/bin/sh

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./qv.sh asm.meryl read.meryl peak threads"
	echo
	echo "  asm.meryl     meryl db of the assembly"
	echo "  read.meryl    meryl db of the read set"
	echo "  peak          haploid peak to normalize kmer multiplicity of the read.meryl"
	echo "  threads       number of threads to use"
	echo
	echo "Output (stdout) format: asm.meryl <tab> excessive <tab> total <tab> qv <tab> error"
	echo "    asm.meryl   name of the provided meryl db of the assembly"
	echo "    excessive   excessive number of kmers found in the assembly"
	echo "    total       total number of kmers in the assembly"
	echo "    qv          estimated base level quality value"
	echo "    error       error rate used to calculate the qv"
	echo
	echo "Aug. 7, 2020"
	echo "Giulio Formenti gformenti@rockefeller.edu"
	echo "Arang  Rhie     arrhie@gmail.com"
	echo
	exit -1;
fi

asm=$1
red=$2
p=$3
threads=$4
if [ -z $threads ]; then
  threads=8
fi

meryl threads=$threads intersect $red $asm output ${red}_and_$asm
meryl threads=$threads divide-round $p ${red}_and_$asm output kr.meryl
meryl threads=$threads subtract $asm kr.meryl output ka-extra.meryl

k=`meryl threads=$threads print $asm | head -n 2 | tail -n 1 | awk '{print length($1)}'`

# echo "# QV statistics for $asm"
ASM_ONLY=`meryl statistics ka-extra.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
TOTAL=`meryl statistics $asm  | head -n4 | tail -n1 | awk '{print $2}'`
ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR"
