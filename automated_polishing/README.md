# Automated Polishing

This README contains details about applying the automated polishing on general genome assemblies.
The [original script](https://github.com/arangrhie/T2T-Polish/blob/master/automated_polishing/automated-polishing-legacy.sh) used 
in [McCartney et al, 2021](https://doi.org/10.1101/2021.07.02.450803) has been updated to use the Merfin release version and contains minor corrections to
Winnowmap alignments.

## Dependencies 
* [Winnowmap2](https://github.com/marbl/Winnowmap)
* [Falconc available in pbipa package](https://github.com/bio-nim/pb-falconc/releases)
* [Racon (liftover branch)](https://github.com/isovic/racon/tree/liftover)
* [Meryl v1.3](https://github.com/marbl/meryl)
* [Merfin v1.0](https://github.com/arangrhie/merfin)
* [Samtools](https://github.com/samtools/samtools)
* [BCFtools](https://github.com/samtools/bcftools)

(Optional)
* [Minimap2 and paftools](https://github.com/lh3/minimap2)

## How to run
```
automated-polishing.sh <num_threads> <num_iterations> <in_draft_fasta> <in_reads> <in_readmers> <out_prefix>

  num_threads    Number of threads to use.
  num_iterations Number of polishing iterations to perform.
  in_draft_fasta Path to the input FASTA/FASTQ containing draft sequences for polishing.
  in_reads       Path to the input reads file, in FASTA/FASTQ format.
  in_readmers    Path to a Meryl database of Illumina mers.
  out_prefix     Prefix of the output files. A folder will automatically be created if it does not exist.
```

## Description

This script automatically launches Winnowmap2 to align HiFi reads, falconc, Racon, and Merfin to call and filter polishing edits,
and runs bcftools to generate the polished consensus in one iteration.
Optionally, edits compared to the initial pre-polished consensus is driven by assembly-to-assembly alignments using minimap2 and paftools.

### 1. HiFi read alignment with Winnowmap2

Prior to polishing automated-polishing.sh initially aligns input HiFi reads to the draft assembly using Winnowmap2,
a tool developed to produce highly accurate repeat-aware alignments, using `--MD -W bad_mers.txt -ax map-pb`.

```
# Collect repetitive 15-mers
meryl count k=15 ${in_draft} output merylDB
meryl print greater-than distinct=0.9998 merylDB > ${out_winnowmap_bam}.repetitive_k15.txt

# Align
winnowmap --MD -W ${out_winnowmap_bam}.repetitive_k15.txt -ax map-pb -t ${threads} ${in_draft} ${in_dataset} | \
  samtools view -Sb > ${out_winnowmap_bam}

# Sort
samtools sort --threads ${threads} -o ${out_winnowmap_sorted_bam} ${out_winnowmap_bam}
```

The unpolished draft assembly is used as the target in the first iteration, while every following iteration uses the polished output of the previous stage
as the input target.

### 2. Filtering bam with Falconc

The output alignments are filtered to remove all secondary alignments and alignments with excessive clipping,
using `falconc bam-filter-clipped` tool in the `pbipa` Bioconda package using options `falconc bam-filter-clipped -t -F 0x104`
where `-F` is the filter flag. By default, maximum clipping on either left or right side of an alignment is set to 100bp,
but it is only applied if the alignment is located at least 25bp from a target sequence end
(to prevent clipping due to contig which could otherwise cause false alignment filtering).

```
falconc bam-filter-clipped -t -F 0x104 \
  --input-fn ${out_winnowmap_sorted_bam} \
  --output-fn ${out_falconc_sam} \
  --output-count-fn ${out_falconc_sam}.filtered_aln_count.txt
```

### 3. Collect polishing edits with Racon

The filtered alignments produced are then used as input to Racon, here the Racon `liftover` branch is utilised.
This is an extension of the `master` branch of Racon with two custom features:
* BED selection of regions for polishing 
* logging the changes introduced to the draft sequences to produce the polished output (in VCF, PAF or optionally SAM format)

Racon is run with default options except for two new logging options `-L out_prefix` and `-S`,
which store the liftover information between the input and output sequences.
We only need the out_prefix.vcf, which is generated with

```
racon -t ${threads} \
  ${in_dataset} \
  ${out_falconc_sam} \
  ${in_draft} \
  -L ${out_racon_fasta} \
  > ${out_racon_fasta}
```

### 4. Filter edits with Merfin

Merfin is run using `-polish` function and `-peak 160.7` which should be adjusted by the user to accommodate their genome
and is calculated from kmer histogram or computed using tools such as [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0).
Refer to [Merfin](https://github.com/arangrhie/merfin) and the [WiKi page](https://github.com/arangrhie/merfin/wiki/Best-practices-for-Merfin)
for more details on how to obtain `-peak kcov` and `-prob lookup_table.txt` using the `--fitted_hist` option in GenomeScope2 and the read kmer database.

```
merfin -polish \
  -sequence ${in_draft} \
  -seqmers ${out_meryl} \
  -readmers ${in_readmers} \
  -peak 106.7 \
  -vcf ${out_racon_fasta}.vcf \
  -output ${out_merfin} \
  -threads ${threads}
```

### 5. Generate consensus with bcftools

Finally, consensus is called using bcftools and to generate the vcf of variants changed from the draft assembly to the polished assembly.
```
bcftools view -Oz ${out_merfin}.polish.vcf > ${out_merfin}.polish.vcf.gz
bcftools index ${out_merfin}.polish.vcf.gz
bcftools consensus ${out_merfin}.polish.vcf.gz -f ${in_draft} -H 1 > ${out_consensus}
```

### 6. Difference between the initial vs. final assembly (Optional)
Additionally, for benchmark purposes, the difference between the initial draft assembly vs. the final polished assembly could be obtained
using minimap2 using the `--cs asm5` function, sorted, and input to paftools for vcf generation.

```
minimap2 -cx asm5 --cs ${in_draft} ${out_consensus} > ${out_consensus}.minimap2.paf
sort -k6,6 -k8,8n ${out_consensus}.minimap2.paf > ${out_consensus}.minimap2.sorted.paf
paftools call -f ${in_draft} ${out_consensus}.minimap2.sorted.paf > ${out_consensus}.paftools.vcf
```

