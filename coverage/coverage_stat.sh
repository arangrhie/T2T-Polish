#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo
  echo "Usage: ./mean_sd.sh in.wig [exclude.bed]"
  echo
  echo "Obtain coverage statistics"
  echo "  in.wig       output from $tools/T2T-Polish/paf_util/pafToCovWig.jar"
  echo "  exclude.bed  region to exclude when calculating coverage stats"
  echo
  exit -1
fi

wig=$1
exc=$2
out=cov

if [[ -f $exc ]]; then
  module load bedtools
  java -jar -Xmx1g $tools/T2T-Polish/paf_util/wigToBed.jar $wig | bedtools subtract -A -a - -b $exc > $out.bed
else
  java -jar -Xmx1g $tools/T2T-Polish/paf_util/wigToBed.jar $wig > $out.bed
fi

cat $out.bed | java -jar -Xmx1g $tools/T2T-Polish/paf_util/txtColumnSummary.jar 4 - | sort -n -k1 - > $out.hist

cat $out.hist | awk '{sum += $1*$2; total+=$2} END {printf "%.2f\n", sum/total}' > $out.mean.txt
COVERAGE_MEAN=`cat $out.mean.txt`
echo "Mean: $COVERAGE_MEAN"

#  Calculate SD for coverage < 2.5 x Mean
cat $out.hist | awk -v COVERAGE_MEAN=$COVERAGE_MEAN '$1 < COVERAGE_MEAN*2.5 {sum += ((($1 - COVERAGE_MEAN)^2)*$2); total+=$2} END {printf "%.2f\n", sqrt(sum/total)}' > $out.sd_adj.txt
SD=`cat $out.sd_adj.txt`
echo "SD Adjusted from cov < 2.5 x Mean: $SD"

# Output high / low cutoffs
echo "$COVERAGE_MEAN $SD" | awk '{print $1*2.5}' > high_cutoff.txt
echo "$COVERAGE_MEAN $SD" | awk '{print $1/4}' > low_cutoff.txt
high=`cat high_cutoff.txt`
low=`cat low_cutoff.txt`
echo "Low cutoff: $low"
echo "High cutoff: $high"
echo

#  Calculate SD
cat $out.hist | awk -v COVERAGE_MEAN=$COVERAGE_MEAN '{sum += ((($1 - COVERAGE_MEAN)^2)*$2); total+=$2} END {printf "%.2f\n", sqrt(sum/total)}' > $out.sd.txt

SD=`cat $out.sd.txt`
echo "SD: $SD"

#  MAD = median(|X - X^|)
total=`cat $out.hist | awk '{sum+=$2} END {print sum}'`
#echo "total=$total"
med=`echo $total | awk '{printf "%.0f\n", $total/2}'`
#echo "med=$med"

cat $out.hist | awk -v med=$med -v isDone=0 '{idx+=$2; if (isDone>0) {next;} else if (idx > med) {print $1; isDone=1;}}' > $out.med.txt
MED=`cat $out.med.txt`
echo "MED: $MED"

cat $out.hist | awk -v MED=$MED '{dev=$1-MED; if (dev<0) dev*=-1; for (i=0; i<$2; i++) { print dev;}}' | sort -n -k1 | sed -n ${med}p - > $out.mad.txt
MAD=`cat $out.mad.txt`
echo "MAD: $MAD"

THET=`echo "$MAD" | awk '{print 1.4826*$1}'`
echo $THET > $out.theta.txt

echo "Theta: $THET"
