#!/bin/bash

if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./make_rDNA_patch.sh species gap-size"
  echo "  species  species id. REQUIRED"
  echo "  gap-size gap size. 1mb or 5kb. default=1mb. OPTIONAL"
  exit -1
fi

module load pangenome
module load seqtk
module load ucsc
module load bedtools

set -x
ver=v1.4.1r
sp=$1
gs=$2
if [[ -z $2 || "$2" == "?" ]]; then
  gs="1mb"
fi

if [[ -z $1 ]]; then
  echo "no species found. exit."
  exit -1
fi

sp_ver=${sp}_$ver

## Create patch vcf
echo "##fileformat=VCFv4.2" > $sp_ver.header.vcf
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> $sp_ver.header.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample" >> $sp_ver.header.vcf

# Check for target.list
if [[ ! -s target.list ]]; then
echo "No target.list found."
exit -1
fi

# Per-coordinate region
N=`wc -l target.list | awk '{print $1}'`

for i in 3 # $(seq 1 $N)
do
  chr_hap=`sed -n ${i}p target.list | awk '{print $1}'`
  target=`sed -n ${i}p target.list | awk '{print $2}'`
  telo=`sed -n ${i}p target.list | awk '{print $3}'`

  if [[ $telo == "clip" ]]; then
    target1=`echo $target | awk -F "|" '{print $1}'`
    samtools view ${chr_hap}_fix_to_$ver.bam |\
      java -jar -Xmx1g $tools/T2T-Polish/patch/samRefPos2QryPos.jar - $target1 > $chr_hap.patch.bed
    target2=`echo $target | awk -F "|" '{print $2}'`
    samtools view ${chr_hap}_fix_to_$ver.bam |\
      java -jar -Xmx1g $tools/T2T-Polish/patch/samRefPos2QryPos.jar - $target2 >> $chr_hap.patch.bed
  else
    samtools view ${chr_hap}_fix_to_$ver.bam |\
    java -jar -Xmx1g $tools/T2T-Polish/patch/samRefPos2QryPos.jar - $target > $chr_hap.patch.bed
  fi

  cat $chr_hap.patch.bed
  rRegion=`cat $chr_hap.patch.bed |\
    awk '{ if (NR==1) { rSeq=$1; rStart=$2; } \
           if (NR==2) { rEnd=$3 } \
         } END {print rSeq":"(rStart+1)"-"rEnd}'`
  ## modify if we are including telomere
  if [[ $telo == "telo" ]]; then
    rRegion=`echo $target | sed 's/,//g'`
  fi
  echo $rRegion

  samtools faidx $sp_ver.$chr_hap.fa $rRegion |\
    awk 'NR!=1' | tr -d '\n' > ref.seq

  ## left side of the gap
  samtools faidx $chr_hap.1.fa
  LEN=`cat $chr_hap.1.fa.fai | awk '{print $2}'`

  ## Take the 1st line. do not extend with "clip"
  if [[ $telo == "clip" ]]; then
  echo "clip to exact matches"
  qRegion1=`sed -n 1p $chr_hap.patch.bed |\
    awk -v LEN=$LEN '{\
      qSeq=$4; qStart=$5; qEnd=$6; \
      if ($NF=="-") { qEnd=(LEN - $5); qStart=(LEN-$6); } \
      } END { print qSeq":"(qStart+1)"-"qEnd }'`
  else
  qRegion1=`sed -n 1p $chr_hap.patch.bed |\
    awk -v LEN=$LEN '{\
     qSeq=$4; qStart=$5; qEnd=LEN; \
     if ($NF=="-") { qEnd=(LEN - qStart); qStart=0; } \
     } END { print qSeq":"(qStart+1)"-"qEnd }'`
  fi
  echo $qRegion1

  samtools faidx $chr_hap.1.fa $qRegion1 > qry.1.fa
  strand=`sed -n 1p $chr_hap.patch.bed | awk '{print $NF}'`
  if [[ "$strand" == "-" ]]; then
  echo "- : rc qry.1.fa"
  seqtk seq -r qry.1.fa > qry.1.rv.fa
  mv qry.1.rv.fa qry.1.fa
  fi

  ## right side of the gap
  samtools faidx $chr_hap.2.fa
  LEN=`cat $chr_hap.2.fa.fai | awk '{print $2}'`

  if [[ $telo == "clip" ]]; then
    echo "clip to exact matches"
    qRegion2=`sed -n 2p $chr_hap.patch.bed |\
           awk -v LEN=$LEN '{ qSeq=$4; qStart=$5; qEnd=$6; \
             if ($NF=="-") { qEnd=(LEN - $5); qStart=(LEN-$6); } \
           } END { print qSeq":"(qStart+1)"-"qEnd }'`
  else
    qRegion2=`sed -n 2p $chr_hap.patch.bed |\
         awk -v LEN=$LEN '{ qSeq=$4; qStart=0; qEnd=$6; \
           if ($NF=="-") { qStart=(LEN - qEnd); qEnd=LEN; } \
         } END { print qSeq":"(qStart+1)"-"qEnd }'`
  fi
  echo $qRegion2
  samtools faidx $chr_hap.2.fa $qRegion2 > qry.2.fa
  strand=`sed -n 2p $chr_hap.patch.bed | awk '{print $NF}'`
  if [[ "$strand" == "-" ]]; then
    echo "- : rc qry.2.fa"
    seqtk seq -r qry.2.fa > qry.2.rv.fa
    mv qry.2.rv.fa qry.2.fa
  fi

  ## Append gap
  cat qry.1.fa | awk 'NR!=1' | tr -d '\n' > qry.1.seq
  cat qry.2.fa | awk 'NR!=1' | tr -d '\n' > qry.2.seq
  cat qry.1.seq ../gap.$gs.seq qry.2.seq > qry.seq

  set +x
  chr=`head -n1 $chr_hap.patch.bed | awk '{print $1}'`
  pos=`head -n1 $chr_hap.patch.bed | awk '{print ($2+1)}'`
  ref=`cat ref.seq`
  alt=`cat qry.seq`

  ## Per $chr_hap vcf
  cat $sp_ver.header.vcf > $sp_ver.$chr_hap.rDNA_patch.vcf
  echo -e "$chr\t$pos\t.\t$ref\t$alt\t1\t.\t.\tGT\t1/1" >> $sp_ver.$chr_hap.rDNA_patch.vcf
  set -x

  ## At the end, gz and make consensus
  bcftools view --no-version -Oz --threads $SLURM_CPUS_PER_TASK $sp_ver.$chr_hap.rDNA_patch.vcf > $sp_ver.$chr_hap.rDNA_patch.vcf.gz
  bcftools index $sp_ver.$chr_hap.rDNA_patch.vcf.gz

  ## Test consensus
  bcftools consensus -c $sp_ver.$chr_hap.rDNA_patched.chain -f $sp_ver.$chr_hap.fa -HA $sp_ver.$chr_hap.rDNA_patch.vcf.gz > $chr_hap.rDNA_patched.fa
  samtools faidx $chr_hap.rDNA_patched.fa

  $MERQURY/eval/qv.sh ../../meryl/${sp}_Hybrid.k31.meryl $chr_hap.rDNA_patched.fa $chr_hap.rDNA_patched

  ## Test by aligning patched sequence to $sp_ver.$chr_hap.fa - doesn't work, sequence too diverged
  echo ">$chr_hap.rDNA_patched.qry" > $chr_hap.rDNA_patched.qry.fa
  cat qry.seq >> $chr_hap.rDNA_patched.qry.fa
  samtools faidx $chr_hap.rDNA_patched.qry.fa

  ## Append 5kb to the target, and liftOver to the new coords
  echo $rRegion | sed 's/,//g' | awk -F ":" '{print $1"\t"$2}' |\
    awk -F "-" '{print $1"\t"$NF}' | \
    awk '{print $1"\t"($2-1-500000)"\t"($3+500000)}' > tmp_target.bed
  liftOver -ends=500 tmp_target.bed $sp_ver.$chr_hap.rDNA_patched.chain \
    tmp_lifted.bed tmp_failed.txt

  if [[ -s tmp_lifted.bed ]]; then
    echo "$rRegion successfully lifted to $lRegion"
    lRegion=`cat tmp_lifted.bed | awk '{print $1":"$2+1"-"}'` # Lifted region
    else
    echo "$rRegion failed to lift over. Taking arbitrary +1Mb for gap"
    lRegion=`cat tmp_target.bed | awk '{print $1":"($2-1)"-"$3+1000000}'`
  fi

  samtools faidx $chr_hap.rDNA_patched.fa $lRegion \
    > $chr_hap.rDNA_patched.qry.fa
  samtools faidx $chr_hap.rDNA_patched.qry.fa

  out=${chr_hap}_rDNA_patched_to_${ver}
  echo $out

  wfmash --no-split -ad -t$SLURM_CPUS_PER_TASK ${sp_ver}.$chr_hap.fa \
    $chr_hap.rDNA_patched.qry.fa -s50000 -p95 > $out.sam
  # $chr_hap.rDNA_patched.fa -s50000 -p95 > $out.sam
  samtools sort -@$SLURM_CPUS_PER_TASK -T $out.tmp -O bam \
    $out.sam > $out.bam
  samtools index $out.bam
  echo

done

samtools merge -@$SLURM_CPUS_PER_TASK -f -o rDNA_patched_to_${ver}.bam chr*_rDNA_patched_to_${ver}.bam
samtools index rDNA_patched_to_${ver}.bam

bcftools concat -a -D --threads $SLURM_CPUS_PER_TASK \
  --no-version -Oz -o $sp_ver.rDNA_patch.vcf.gz \
  $sp_ver.chr*.rDNA_patch.vcf.gz
bcftools index $sp_ver.rDNA_patch.vcf.gz
