#!/bin/bash

# ED_Fig. b-d (left): Merqury spectra-cn plots

# $FASTA: v0.9 fasta.
# Available at https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v0.9.fasta

FASTA=chm13v0.9.fasta

# Run Merqury for each read dbs: HiFi, Illumina, and Hybrid
#  Illumina.k21.meryl: Download and extract https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/qc/IlluminaPCRfree.k21.meryl.tar.gz
#  HiFi.k21.meryl: Download and extract https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/qc/hifi20kb.k21.meryl.tar.gz
#  Hybrid.k21.meryl: Download and extract https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/qc/hybrid.k21.meryl.tar.gz

for PLATFORM in HiFi Illumina Hybrid
do
  mkdir -p $PLATFORM
  cd $PLATFORM
  $MERQURY/merqury.sh ../$PLATFORM.k21.meryl ../$FASTA $PLATFORM
  cd ../
done
# This generates histograms shown on ED Fig.1 b-d (Left) and chm13v0.9_only.bed under each $PLATFORM.

# ED_Fig. b-d (right): missing k-mers per chromosome
# $PATCHES: Region patched with ONT-version of assembly
# Available at https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/qc/chm13v0.9_patches_only.bed
PATCHES=chm13v0.9_patches_only.bed
MERQURY_ONLY=chm13v0.9_only.bed
OUTPUT=ED_Fig1b-d_missing_kmers.tab

echo -e "CHR\tMISSING\tTYPE\tPLATFORM" > $OUTPUT

for PLATFORM in HiFi Illumina Hybrid
do
  bedtools intersect -u -a $PLATFORM/$MERQURY_ONLY -b $PATCHES  | \
  java -jar -Xmx1g $tools/T2T-Polish/paf_util/txtColumnSummary.jar 1 - |\
  sed 's/chr//g' | awk -v PLATFORM=$PLATFORM '{print $0"\tPatches\t"PLATFORM}' |\
  sort -k1 -n >> $OUTPUT
  bedtools subtract -a $PLATFORM/$MERQURY_ONLY -b $PATCHES      |\
  java -jar -Xmx1g $tools/T2T-Polish/paf_util/txtColumnSummary.jar 1 - |\
  sed 's/chr//g' | awk -v PLATFORM=$PLATFORM '{print $0"\tHiFi_Consensus\t"PLATFORM}' |\
  sort -k1 -n >> $OUTPUT
done

Rscript ED_Fig1b-d.R