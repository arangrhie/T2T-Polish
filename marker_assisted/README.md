# Marker assisted read filtering

## Requires
* meryl-1.0
* python 2.7.3
* java
* samtools
* IGVtools (for generating .tdf coverage file)

## Prepare marker db

To generate one, generate a k-mer counting database with `meryl count` and filter to have the matching “marker” criteria.

Depending on the purpose, this marker set can be generically defined.

As the polishing team aims to generate the most conservative alignment set to validate structural correctness, we define markers as single-copy k-mers, occurring in the expected single-copy range from the PCR-free Illumina read set and unique in the assembly.

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

Markers were generated from Illumina PCR-Free WGS 21-mer counts (IlluminaPCRfree.meryl) and 21-mers of an assembly (20200727.meryl) using the following command lines:
```
# Get k-mers with multiplicity > 42, filter out erroneous k-mers
meryl greater-than 42 IlluminaPCRfree.meryl output IlluminaPCRfree.gt42.meryl

# Get k-mers with multiplicity < 133, filter out > 2.5 copy k-mers (multiplicity at first peak + 1/4 * (second - first peak))
meryl less-than 133 IlluminaPCRfree.gt42.meryl output IlluminaPCRfree.gt42.lt133.meryl
# IlluminaPCRfree.gt42.lt133.meryl is available as IlluminaPCRfree.single.meryl

# Get unique k-mers in the assembly
meryl equal-to 1 20200727.meryl output 20200727.single.meryl

# Intersect to guarantee the single-copy k-mers in the reads are globally unique in the assembly
meryl intersect IlluminaPCRfree.single.meryl 20200727.single.meryl output IlluminaPCRfree.single.20200727.meryl
```

Our markers associated with a released t2t-chm13 assembly will be available on Globus as:

`team-curation/marker_assisted/IlluminaPCRfree.single.YYYYMMDD.meryl.tar.gz`

## Run

This script is composed of three steps:
1.	init.sh: Extract mappings to a `target` sequence, sort by read id, split alignments per 10k read chunks
2.	convert.sh: Sort each read by SAM flag and convert each chunk to a fasta file preserving mapping position and orientation. A job array is submitted with each processing a chunk.
3.	merge.sh: Count markers in each alignment, preserve alignments having the most unique markers when a subsequence is mapped to multiple positions (markers.cram), filter for length by `length_filt` and alignment identity >75% (markersandlength.cram), and generate the .tdf coverage tracks

Submit example:
```
./_submit.sh input.bam chr20 t2t-chm13.20200727.fasta IlluminaPCRfree.single.20200727.meryl 25
```
* input.bam: Any alignment file. We use winnowmap alignments. Sort and index.
* target: Sequence id to generate the alignments. ex. chr20
* assembly.fa: Full assembly set used to generate the input.bam. ex. t2t-chm13.20200602.fasta
* marker.meryl: marker db.
* length_filter: Filter out alignments shorter than this length. in kbp.
