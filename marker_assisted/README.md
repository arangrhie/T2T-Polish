# Marker assisted read filtering

## Requires
* [meryl-1.3](https://github.com/marbl/meryl)
* python 2.7.3 or higher
* java
* samtools

## Prepare marker db

To generate one, generate a k-mer counting database with `meryl count` and filter to have the matching “marker” criteria.

Depending on the purpose, this marker set can be generically defined.

As the goal is to generate the most conservative alignment set to validate structural correctness, we define markers as single-copy k-mers,
occurring in the expected single-copy range from the PCR-free Illumina read set and unique in the assembly.

Note the "single-copy range" is different from the conventional definition of single-copy k-mer multiplicity found in a genome.
Here, we refer to the k-mers that are expected to be found 'globaly once' in a pseudo-haploid assembly as "single-copy",
where k-mers from the heterozygous (1-copy in the genome) and homozygous (2-copy in the genome) regions are expected to be found once.

As CHM13 haplotypes are nearly identical, we found almost no strong signal for the heterozygous 1-copy k-mer peak.
This makes the first visible peak observed in the Illumina k-mer set become the 2-copy peak.
The second peak was nearly doubling the multiplicity found from what was found in the first peak, making it the 4-copy peak.

When evaluating read alignments at heterozygous regions (we still found a few heterozygous sites), 
we confirmed the heterozygous regions matched the expected haploid multiplicity (half the first peak multiplicity), 
and the average read depth (usually called haploid genome coverage) matched the multiplicity of the first peak.

Using this information, we set the upper boundary of the single-copy k-mer multiplicity as the 2.5-copy multiplicity, which can be inferred from the first (2-copy) and second (4-copy) peaks.

We chose k=21 following [Fofanov et al, Bioinformatics, 2004](https://doi.org/10.1093/bioinformatics/bth266) to allow a maximum k-mer collision rate of 0.005, which is close to the Illumina sequencing error rate, from a 3.2Gb human genome size using [this script](https://github.com/marbl/merqury/blob/master/best_k.sh).

Markers were generated from [Illumina PCR-Free WGS](https://github.com/marbl/CHM13#illumina-pcrfree-data) 21-mer counts [chm13.k21.meryl](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/MERFIN_2021/chm13/evaluation/chm13.k21.meryl.tar.gz).

In brief,
```
# Get k-mers with multiplicity > 42, filter out erroneous k-mers
meryl greater-than 42 chm13.k21.meryl output chm13.gt42.meryl

# Get k-mers with multiplicity < 133, filter out > 2.5 copy k-mers (multiplicity at first peak + 1/4 * (second - first peak))
meryl less-than 133 chm13.gt42.meryl output IlluminaPCRfree.single.k21.meryl
```

Unique 21-mers of an assembly can be generated using the following command lines:
```
# Get unique k-mers in the assembly
meryl count k=21 $asm.fasta output $asm.meryl
meryl equal-to 1 $asm.meryl output $asm.1.meryl
```

The final marker set is generated through an intersection:
```
# Intersect to guarantee the single-copy k-mers in the reads are globally unique in the assembly
meryl intersect IlluminaPCRfree.single.k21.meryl $asm.1.meryl output $asm.single.k21.meryl
```

### Markers are available below:
* [IlluminaPCRfree.single.k21.meryl](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/IlluminaPCRfree.single.k21.meryl.tar.gz)
* [IlluminaPCRfree.single.k31.meryl](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/IlluminaPCRfree.single.k31.meryl.tar.gz)
* [IlluminaPCRfree.single.k51.meryl](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/IlluminaPCRfree.single.k51.meryl.tar.gz)

### Markers associated with a released T2T-CHM13 assembly is available below:
#### v0.9
* [chm13v0.9.single.k21.meryl.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/chm13v0.9.single.k21.meryl.tar.gz)

#### v1.0
* [chm13v1.0.single.k21.meryl.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/chm13v1.0.single.k21.meryl.tar.gz)
* [chm13v1.0.single.k51.meryl.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/chm13v1.0.single.k51.meryl.tar.gz)

#### v1.1
* [chm13v1.1.single.k21.meryl.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/chm13v1.1.single.k21.meryl.tar.gz)
* [chm13v1.1.single.k31.meryl.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/chm13v1.1.single.k31.meryl.tar.gz)
* [chm13v1.1.single.k51.meryl.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/marker/chm13v1.1.single.k51.meryl.tar.gz)

## Run

The marker assisted alignment runs per chromosome. Make sure the input alignment bam file is sorted and indexed.

### 1. Filter alignment per chromosome
```
for i in $(seq 1 22) X
do
  filter_by_marker_nosplit.sh $alignment.bam chr$i $asm.fasta single.meryl len_filt
done
```

A brief help message on how to run this script:
```
./filter_by_marker_nosplit.sh
Usage: filter_by_marker.sh alignment target asm marker.meryl len_filt

Filters alignment based on single-copy kmers
  alignment: input bam or cram file
  target: target region (chr) to process
  asm: entire fasta file
  marker.meryl: meryl db containing marker kmers. This scripts generates chr specific markers.
  len_filt: alignment length filter in kb. ex. 1 for 1kb
```

`len_filt` are alignment length filter in kb, to remove spurious short alignments. Below are possible length filters for long-reads (used in [McCartney et al, 2021](https://doi.org/10.1101/2021.07.02.450803)) and for short-reads (used in [Hoyet et al, 2021](https://doi.org/10.1101/2021.07.12.451456)).
* PacBio HiFi: `10`
* PacBio CLR: `1`
* ONT: `25`
* Illumina: `0`

This step counts markers in each alignment, preserves alignments having the most unique markers when a subsequence is mapped to multiple positions (markers.bam), filter for length by `length_filt` and alignment identity >75% (markersandlength.bam).

### 2. Aggregate per-chromosome filtered alignments
```
aggregate.sh $out_prefix
```
This step merges all `\*.markers.bam` and `\*.markersandlength.bam` separately, converts to `.paf`, and generates a coverage-depth `.wig` file.
