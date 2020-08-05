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

As the polishing team aims to generate the most conservative alignment set to validate structural correctness, we define markers as single-copy k-mers, occurring in the expected single-copy range from the read set and unique in the assembly.

Markers were generated from Illumina PCR-Free WGS 21-mer counts (IlluminaPCRfree.meryl) and 21-mers of an assembly (20200727.meryl) using the following command lines:
```
# Get k-mers with multiplicity > 42, filter out erroneous k-mers
meryl greater-than 42 IlluminaPCRfree.meryl output IlluminaPCRfree.gt42.meryl

# Get k-mers with multiplicity < 133, filter out >1.5 copy k-mers
meryl less-than 133 IlluminaPCRfree.gt42.meryl output IlluminaPCRfree.gt49.lt133.meryl
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
