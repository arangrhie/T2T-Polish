#!/bin/bash

if [[ "$#" -lt 2 ]]; then
  echo "Usage: ./collect_summary.sh in.paf name"
  exit -1
fi

pri=$1
pri=${pri/.paf/}
name=$2

collect_stat() {
  in_paf=$1
  out_prefix=$2
  WINDOW=$3 # 10000
  platform=$4
  
  perstrand=$5

  if [[ "$perstrand" -eq 1 ]]; then
# Read length
    awk '$5=="+"' $in_paf |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform Len (+)" $WINDOW 2 $out_prefix.readLen.p

    awk '$5=="-"' $in_paf |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform Len (-)" $WINDOW 2 $out_prefix.readLen.n

# MQ
    awk '$5=="+"' $in_paf |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform MQ (+)" $WINDOW 12 $out_prefix.mq.p

    awk '$5=="-"' $in_paf |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform MQ (-)" $WINDOW 12 $out_prefix.mq.n

# Idy
    cut -f1-11 $in_paf |\
    awk '$5=="+" {print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform Idy (+)" $WINDOW 12 $out_prefix.idy.p

    cut -f1-11 $in_paf |\
    awk '$5=="-" {print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform Idy (-)" $WINDOW 12 $out_prefix.idy.n

  else
    cat $in_paf |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform Len" $WINDOW 2 $out_prefix.len

    cat $in_paf |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform MQ" $WINDOW 12 $out_prefix.mq
  
    cut -f1-11 $in_paf |\
    awk '{print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar - "$platform Idy" $WINDOW 12 $out_prefix.idy

# Strand
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToStrandedWig.jar $in_paf "$platform +/all" $WINDOW norm > $out_prefix.strand.wig
  fi

}

collect_clipped() {
  in_paf=$1
  out_prefix=$2
  WINDOW=$3
  platform=$4

  perstrand=$5

# Clipped bases > 100bp
  java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar $in_paf "$platform Clipped (%)" $WINDOW norm 100 > $out_prefix.clip_norm.wig
  java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar $in_paf "$platform Clipped" $WINDOW abs 100 > $out_prefix.clip_abs.wig

# Clipped per strand
  if [[ $perstrand -eq 1 ]]; then
    cut -f1-11 $in_paf |\
    awk '$5=="+" {print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar - "$platform Clipped + (%)" $WINDOW norm 100 > $out_prefix.p.clip_norm.wig

    cut -f1-11 $in_paf |\
    awk '$5=="+" {print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar - "$platform Clipped +" $WINDOW abs 100 > $out_prefix.p.clip_abs.wig

    cut -f1-11 $in_paf |\
    awk '$5=="-" {print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar - "$platform Clipped - (%)" $WINDOW norm 100 > $out_prefix.n.clip_norm.wig

    cut -f1-11 $in_paf |\
    awk '$5=="-" {print $0"\t"100*$10/$11}' - |\
    java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar - "$platform Clipped -" $WINDOW abs 100 > $out_prefix.n.clip_abs.wig
  fi
}

collect_coverage() {
  in_paf=$1
  out_prefix=$2
  WINDOW=$3
  platform=$4
  perstrand=$5

  if [[ $perstrand -eq 1 ]]; then
    cut -f1-11 $in_paf |\
    awk '$5=="+"' |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToCovWig.jar - "$platform (+)" $WINDOW > $out_prefix.cov.p.wig

    cut -f1-11 $in_paf |\
    awk '$5=="-"' |\
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToCovWig.jar - "$platform (-)" $WINDOW > $out_prefix.cov.n.wig
  else
    java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToCovWig.jar $in_paf "$platform" $WINDOW > $out_prefix.cov.wig
  fi
}


# Collect read-length, mq, idy tracks
collect_stat $pri.paf $pri 10000

# Per strand stat tracks
collect_stat $pri.paf $pri 10000 "$name" 1

# Clipped
collect_clipped $pri.paf $pri.w1k 1024 "$name"

# Coverage
collect_coverage $pri.paf $pri 1024 "$name"

# Per strand coverage
collect_coverage $pri.paf $pri 1024 "$name" 1




