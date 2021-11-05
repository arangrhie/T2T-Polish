#!/bin/bash

if [[ "$#" -lt 1 ]]; then
	echo "Usage: ./aggregate.sh out-prefix"
	exit -1
fi

cpus=$SLURM_CPUS_PER_TASK
out=$1

aggregate() {
	prefix=$1
	out=$2
	ls chr*.$prefix.bam > $prefix.list

	echo "
	Samtools Merge"
	samtools merge -b $prefix.list -O bam -@$cpus $out.$prefix.bam
	samtools index $out.$prefix.bam

	echo "
	Convert to paf"
	$tools/T2T-Polish/coverage/sam2paf.sh $out.$prefix.bam $out.$prefix.paf
	java -jar -Xmx4g $tools/T2T-Polish/paf_util/pafToCovWig.jar $out.$prefix.paf "$out" 1024 > $out.$prefix.cov.wig
  pigz $out.$prefix.paf
}

aggregate markersandlength $out
aggregate markers $out

