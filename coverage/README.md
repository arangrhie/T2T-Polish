# Coverage analysis

Once the alignments are ready in paf format, various tracks and statistics can be generated.

## Dependencies
* [seqtk](https://github.com/lh3/seqtk)
* [seqrequester](https://github.com/marbl/seqrequester)
* [asset](https://github.com/dfguan/asset/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* java

## Prerequisites
All scripts requires `$tools` variable to be set prior to running.  
The `$tools` is the dir path where this github repository has been cloned under.  
`seqrequester` and `asset` should be also installed under `$tools`.

## Quick start
```
# This step needs to be run only once per assembly.
mkdir pattern && cd pattern
$tools/T2T-Polish/coverage/init.sh /path/to/asm_ver.fasta
ln -s /path/to/merqury/asm_ver_only.bed asm_ver.error.bed  # symlink asm_only.bed file from Merqury. Hybrid kmers recommended.

cd /path/to/mappings

# Run this step on each in.paf file per platform.
$tools/T2T-Polish/coverage/issues.sh in.paf name ver platform pattern

# Load *.cov.wig and *.issues.bed on IGV to visually inspect coverage and related issues.
```

Here is the discription of `init.sh` and `issues.sh`:
```
$tools/T2T-Polish/coverage/init.sh
Usage: ./init.sh <asm.fasta>
  <asm.fasta>: absolute path to asm.fasta,
               expecting to have asm.fasta.fai in the same path

$tools/T2T-Polish/coverage/issues.sh
Usage: ./issues.sh in.paf name ver platform pattern
  in.paf   : input paf file
  name     : name to appear on the .wig files
  ver      : assembly version
  platform HiFi, ONT or Hybrid
  pattern  path to pattern folder

Required: pattern/ver/, pattern/ver.bed, pattern/ver.error.bed, pattern/ver.exclude.bed, pattern/ver.telo.bed
```

This runs the following scripts:
* collect_summary.sh
* coverage_stat.sh
* low_support.sh

`low_support.sh` requires additional files. Adjust the script to match the file path. See [prerequisites](#prerequisites-1).

## Collect Summary

```
Usage: ./collect_summary.sh in.paf name
```
By default, using a given paf file, this script generates various summary `.wig` tracks:
* collect_stat: all and per-strand read length, MQ, Idy, and strand ratio in every 10kb window
* collect_coverage: coverage and clipped (>100bp) read tracks in every 1024 bp

Internally, this launches several paf tools that could be run independently.<br>
The track will be displayed with as `name`. Use quotes (`" "`) if the name contains spaces.


### pafToMinAvgMedMaxWig.jar

Generate `minimum`, `average`, `median`, `maximum` `.wig` tracks from values in `col-idx` in every `span` bp. 
```
java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToMinAvgMedMaxWig.jar <in.paf> <name> <span> <col-idx> <out-prefix>
  <name>       : name of this track. String
  <span>       : span of the interval. INT
  <col-idx>    : column to collect stats. DEFAULT=2
  <out-prefix> : output prefix
  output       : <out-prefix.type.wig>. Four .wig formatted files are generated; with type = min avg med and max.
```

### pafToStrandedWig.jar

Generate `+`, `-`, or `norm`alized count tracks by the total reads in every `span` bp. Output `.wig` file is generated in standard out.
```
java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToStrandedWig.jar <in.paf> <name> <span> <type>
  <name> : name of this track. String
  <span> : span of the interval. INT
  <type> : + | - | norm. norm = normalized + counts by total reads
  stdout : .wig format.
```

### pafToCovClippedWig.jar

Generate coverage (read counts), `abs`olute and `norm`alized counts of clipped and total read tracks in every `span` bp; if the soft-clipped or hard-clipped bases are over `min-clipped`.

```
Usage: java -jar -Xmx8g pafToCovClippedWig.jar <in.paf> <name-prefix> <span> <out-prefix> [min-clipped]
  name-prefix name prefix of the tracks. String
  span        span of the interval. INT
  out-prefix  output prefix.
  min-clipped minimum num. of clipped bases. DEFAULT=100. OPTIONAL.
```

## Coverage mean and SD

Using the `.cov.wig` file from `collect_summary.sh`, obtain mean and SD.
```
Usage: ./mean_sd.sh in.wig [exclude.bed]

Obtain coverage statistics
  in.wig       output from /data/Phillippy/tools/T2T-Polish/paf_util/pafToCovWig.jar
  exclude.bed  region to exclude when calculating coverage stats (optional)
```
This script outputs the followings:
* cov.mean.txt   : Coverage mean
* high_cutoff.txt: Mean x 2.5
* low_cutoff.txt : Mean / 4
* cov.out.sd.txt : SD from all coverage
* cov.sd_adj.txt : Adjusted SD using coverage under 2.5 x mean
* cov.med.txt    : Median(|X - X^|)
* cov.mad.txt    : [Median absolute deviation (MAD)](https://en.wikipedia.org/wiki/Median_absolute_deviation)
* cov.theta.txt  : 1.4826 x MAD

`high_cutoff.txt` and `low_cutoff.txt` are used to detect low coverage regions (likely mis-assembly) and high coverage regions (likely collapses). These cutoffs could be adjusted depending on the coverage.

## Low support

This script generates the `issues.bed` file, reported in [CHM13-Issues](https://github.com/marbl/CHM13-issues) repo.

```
Usage: ./low_support.sh in.paf ver platform pattern
  in.paf    sam2paf
  ver       assembly version prefix
  platform  HiFi, ONT or Hybrid
  pattern   path to pattern folder
```
Adjust the file paths below before running this script.

### Prerequisites

`issues.sh` requires several files under `pattern`.  
`asm`, `telo`, `exclude`, `pattern` will be generated with [`init.sh`](init.sh).

* asm    : `ver`.bed; Bed file containing regions as chromosome (scaffold) length.
* telo   : `ver`.telo.bed; Both ends of each chromosome. Will be excluded from flagging as lower coverage is expected due to natural sequencing coverage dropouts. Could be an empty file.
* error  : `ver`.error.bed; Regions of known consensus base errors. [Merqury](https://github.com/arangrhie/T2T-Polish/tree/master/merqury) asm only track.
* exclude: `ver`.exclude.bed; Any additional regions to exclude, e.g. known collapses or gaps. Could be an empty file.
* pattern: pattern/`ver` folder; containing microsatellite patterns.

### Hybrid mode
Asset runs without excluding any bases at the beginning or end of a read alignment while collecting coverage (`-l 0`).  
Clipping threshold is set to `5`, flagging regions when > 5 reads or > 5% of the reads have clippings.

### HiFi mode
Asset runs without excluding any bases at the beginning or end of a read alignment while collecting coverage (`-l 0`).  
Clipping threshold is set to `10`, flagging regions when > 10 reads or > 10% of the reads have clippings.

### ONT mode
Asset runs with default options.  
Clipping threshold is set to `15`, flagging regions when > 15 reads or > 15% of the reads have clippings.

### Output file

Output prefix follows the in.paf file (without the .paf).
* prefix.issues.bed
* prefix.clipped.bed

#### prefix.issues.bed
Regions are reported with one of the following categories:

| Label | Description | R,G,B | Color|
| :--- | :--- | :---: | :---: |
| Low | Low coverage | 204,0,0 | red |
| Low_Qual | Low coverage from lower consensus quality | 204,0,0 | red |
| Low_[GA/TC\|GC\|AT] | Low coverage associated with > 80% microsatellites | 153,153,255 | light purple |
| Clipped | Clipped regions | 153,153,153 | Grey |
| High | High coverage | 153,102,255 | purple |

Precedence takes place in the following order:
* Low_Qual: Low coverage, associated with an error k-mer
* Low_GA/TC|GC|AT: Low coverage, associated with possible sequencing bias. GA|TC > GC > AT.
* Low
* High
* Clipped

#### prefix.clipped.bed
Clipped regions not reported because of the overlapping low or high coverage are reported in this file.
