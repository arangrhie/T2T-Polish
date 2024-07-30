#!/bin/bash

### Need the following files:
# $sp_ver.$chr_hap.fa     target reference sequence
# ${chr_hap}_telo_to_$ver.bam  alignment of the patch to target
# $chr_hap.telo.fa        telomere patch sequence
# target.telo.list        chr_hap chr_hap:start-end(target) telo(p or q)


if [[ "$#" -lt 1 ]]; then
  echo "Usage: ./make_telo_patch.sh species"
  echo "  species  species id. REQUIRED"
  exit -1
fi

module load pangenome
module load seqtk
module load ucsc
module load bedtools

set -x
ver=v1.4.1r
sp=$1

if [[ -z $1 ]]; then
  echo "no species found. exit."
  exit -1
fi

set -x
ver=v1.4.1r
sp_ver=${sp}_$ver

## Create patch vcf
if [[ ! -s  $sp_ver.header.vcf ]]; then
echo "##fileformat=VCFv4.2" > $sp_ver.header.vcf
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> $sp_ver.header.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample" >> $sp_ver.header.vcf
fi

TLIST=target.telo.list
# Check for $TLIST
if [[ ! -s $TLIST ]]; then
echo "No $TLIST found."
exit -1
fi

# Per-coordinate region
N=`wc -l $TLIST | awk '{print $1}'`

for i in $(seq 1 $N)
do
  chr_hap=`sed -n ${i}p $TLIST | awk '{print $1}'`
  target=`sed -n ${i}p $TLIST | awk '{print $2}'`
  telo=`sed -n ${i}p $TLIST | awk '{print $3}'` # p or q

  samtools view ${chr_hap}_telo_to_$ver.bam |\
    java -jar -Xmx1g ~/codes/samRefPos2QryPos.jar - $target > $chr_hap.telo_$telo.bed

  if [[ "$telo" == "p" ]]; then
    cat $chr_hap.telo_$telo.bed
    rRegion=`cat $chr_hap.telo_$telo.bed |\
             awk '{ print $1":1-"$3 }'`
  elif [[ "$telo" == "q" ]]; then
    LEN=`cat $sp_ver.$chr_hap.fa.fai | awk '{print $2}'`
    cat $chr_hap.telo_$telo.bed
    rRegion=`cat $chr_hap.telo_$telo.bed |\
             awk -v LEN=$LEN '{ print $1":"($2+1)"-"LEN }'`
  fi
  echo $rRegion

  samtools faidx $sp_ver.$chr_hap.fa $rRegion | awk 'NR!=1' | tr -d '\n' > ref.seq

# Telomere patch sequence
  samtools faidx $chr_hap.telo.fa
  LEN=`cat $chr_hap.telo.fa.fai | awk '{print $2}'`

## Extend for p and q
  if [[ $telo == "p" ]]; then
  qRegion=`cat $chr_hap.telo_$telo.bed |\
           awk -v LEN=$LEN \
           '{ qSeq=$4; qStart=0; qEnd=$6; \
              if ($NF=="-") { qStart=(LEN-$6); qEnd=LEN; } \
            } END { print qSeq":"(qStart+1)"-"qEnd }'`
  elif [[ $telo == "q" ]]; then
    qRegion=`cat $chr_hap.telo_$telo.bed |\
            awk -v LEN=$LEN \
            '{ qSeq=$4; qStart=$5; qEnd=LEN; \
              if ($NF=="-") {  qStart=0; qEnd=(LEN - $5);} \
             } END { print qSeq":"(qStart+1)"-"qEnd }'`
  fi
  echo $qRegion

  # Reverse complement the fa if it's in - direction
  samtools faidx $chr_hap.telo.fa $qRegion > qry.fa
  strand=`cat $chr_hap.telo_$telo.bed | awk '{print $NF}'`
  if [[ "$strand" == "-" ]]; then
    echo "rc qry.fa"
    seqtk seq -r qry.fa > qry.rv.fa
    mv qry.rv.fa qry.fa
  fi

  # Remove \n and > lines
  cat qry.fa | grep -v '^>' | tr -d '\n' > qry.seq

  set +x
  chr=`cat $chr_hap.telo_$telo.bed | awk '{print $1}'`

  if [[ "$telo" == "p" ]]; then
    pos=1;
  else
    pos=`cat $chr_hap.telo_$telo.bed | awk '{print ($2+1)}'`
  fi
  ref=`cat ref.seq`
  alt=`cat qry.seq`

## Per $chr_hap vcf
  cat $sp_ver.header.vcf > $sp_ver.$chr_hap.telo_patch.vcf
  echo -e "$chr\t$pos\t.\t$ref\t$alt\t1\t.\t.\tGT\t1/1" >> $sp_ver.$chr_hap.telo_patch.vcf
  set -x

## At the end, gz and make consensus
  bcftools view --no-version -Oz --threads $SLURM_CPUS_PER_TASK \
    $sp_ver.$chr_hap.telo_patch.vcf > $sp_ver.$chr_hap.telo_patch.vcf.gz
  bcftools index $sp_ver.$chr_hap.telo_patch.vcf.gz

## Test consensus
  bcftools consensus -c $sp_ver.$chr_hap.telo_patch.chain \
    -f $sp_ver.$chr_hap.fa -HA \
    $sp_ver.$chr_hap.telo_patch.vcf.gz > $chr_hap.telo_patched.fa
  samtools faidx $chr_hap.telo_patched.fa

## Test by aligning $chr_hap.rDNA_test.fa to $sp_ver.$chr_hap.fa
  out=${chr_hap}.${telo}_telo_patched_to_${ver}
  echo $out

## Map patched reference to $ver
  if [[ "$telo" == "p" ]]; then
    samtools faidx $chr_hap.telo_patched.fa $chr_hap:1-800000 > $chr_hap.telo_patched.$telo.fa
  else
    rEnd=`cat $chr_hap.telo_patched.fa.fai | awk '{print $2}'`
    samtools faidx $chr_hap.telo_patched.fa $chr_hap:$(($rEnd-800000))-$rEnd > $chr_hap.telo_patched.$telo.fa
  fi
  samtools faidx $chr_hap.telo_patched.$telo.fa
  wfmash --no-split -ad -t$SLURM_CPUS_PER_TASK ${sp_ver}.$chr_hap.fa \
    $chr_hap.telo_patched.$telo.fa -s50000 -p95 > $out.sam
  samtools sort -@$SLURM_CPUS_PER_TASK -T $out.tmp -O bam \
    $out.sam > $out.bam
  samtools index $out.bam
  echo

#$MERQURY/eval/qv.sh ../../meryl/${sp}_Hybrid.k31.meryl $chr_hap.telo_patched.fa $chr_hap.telo_patched
  $MERQURY/eval/qv.sh ../../meryl/${sp}_Hybrid.k31.meryl $chr_hap.telo_patched.$telo.fa $chr_hap.telo_patched.$telo

done
samtools merge -@$SLURM_CPUS_PER_TASK -f \
  -o telo_patched_to_${ver}.bam chr*_telo_patched_to_${ver}.bam
samtools index telo_patched_to_${ver}.bam

bcftools concat -a -D --threads $SLURM_CPUS_PER_TASK \
  --no-version -Oz -o $sp_ver.telo_patch.vcf.gz \
  $sp_ver.chr*.telo_patch.vcf.gz
bcftools index $sp_ver.telo_patch.vcf.gz

