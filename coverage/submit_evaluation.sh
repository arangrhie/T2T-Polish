#!/bin/sh
#SBATCH --job-name=t2t-eval
#SBATCH --partition=quick
#SBATCH --cpus-per-task=24
#SBATCH --mem=48g
#SBATCH --time=4:00:00
#SBATCH --chdir=.
#SBATCH --output=logs/t2t-eval.%j.log
#SBATCH --error=logs/t2t-eval.%j.log

if [[ "$#" -lt 1 ]]; then
  echo "Usage: $0 asm.dip.fasta.gz"
  echo "Submit post assembly evaluation jobs"
  echo "  asm.dip.fasta.gz: diploid assembly FASTA from T2T-Polish, e.g. assembly_v0.5.dip.fa.gz"
  exit 1
fi

# Input: results folder from T2T-Polish nextflow pipeline with `--run_dv false`
# Output: results/merqury, results/pattern and results/issues

mkdir -p logs

results=$(realpath results)
meryl=$(realpath hybrid.k31.meryl) # hybrid.k31.meryl: meryl database from HiFi and Illumina / Element short-reads, k=31
asm=$(realpath $1)
name=$(basename $asm)
name=`echo $name | sed 's/\.dip.fa.gz$//g'`

echo "== T2T-Polish evaluation for $name== "

module load bedtools/2.31.1
module load ucsc/495

set -e
set -o pipefail
set -x

cd $results
mkdir -p merqury && cd merqury
if [[ -s ${name}_only.bed ]]; then
  echo "Merqury results found."
  cat ${name}.qv
else
  $MERQURY/eval/qv.sh $meryl $asm $name
fi
cd ../

mkdir -p pattern && cd pattern
if [[ -s $name.error.bed && -s $name.bed && -s $name.exclude.bed && -s $name.telo.bed ]]; then
  echo "Pattern files found."
else
  $tools/T2T-Polish/coverage/init.sh $asm $name
  bedtools merge -i ../merqury/${name}_only.bed > ${name}.error.bed
fi
cd ../

mkdir -p issues && cd issues
mkdir -p ont hifi
cd ont
$tools/T2T-Polish/coverage/issues.sh ../../mapping/$name.dip.ont/$name.dip.ont.pri.paf   $name $name ONT ../../pattern
cd ../hifi
$tools/T2T-Polish/coverage/issues.sh ../../mapping/$name.dip.hifi/$name.dip.hifi.pri.paf $name $name HiFi ../../pattern
cd ../

hifi_issues=hifi/$name.dip.hifi.pri.issues.fm.bed
ont_issues=ont/$name.dip.ont.pri.issues.fm.bed

bedtools intersect -u -a $hifi_issues -b $ont_issues > ${name}.issues.bed
bedToBigBed -type=bed9 ${name}.issues.bed ../pattern/${name}.fa.fai ${name}.issues.bb
set -x

echo "Done!"