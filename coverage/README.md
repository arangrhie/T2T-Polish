# Coverage analysis

Once the alignments are ready in paf format, various tracks and statistics can be generated.

## Dependencies
* [Asset](https://github.com/dfguan/asset/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* java

## Prerequisites
All scripts requires `$tools` variable to be set prior to running.
The `$tools` is the dir path where this github repository has been cloned under.

## Quick start
```
Usage: ./issues.sh in.paf name ver platform
  in.paf   : input paf file
  name     : name to appear on the .wig files
  ver      : assembly version
  platform : HiFi or ONT
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
* collect_clipped: clipped reads in every 1024 bp window
* collect_coverage: coverage tracks in every 1024 bp

Internally, this launches several paf tools that could be run independently. The track will be named after `name`. Use quotes (`” “`) if the name contains spaces.


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

### pafToClippedWig.jar

Generate `abs`olute or `norm`alized counts of clipped read tracks in every `span` bp; if the soft-clipped or hard-clipped bases are over `min-clipped`.

```
java -jar -Xmx3g $tools/T2T-Polish/paf_util/pafToClippedWig.jar <in.paf> <name> <span> <type> [min-clipped]
	<name>        : name of this track. String
	<span>        : span of the interval. INT
	<type>        : abs | norm. abs = absolute, norm=normalized by total reads
	[min-clipped] : minimum num. of clipped bases. DEFAULT=100
	stdout: .wig format.
```

### pafToCovWig.jar

Generate coverage tracks (read counts) in every `span` bp.
```
java -jar -Xmx8g $tools/T2T-Polish/paf_util/pafToCovWig.jar <in.paf> <name> <span>
	<name> : name of this track. String
	<span> : span of the interval. INT
	stdout : .wig format.
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
Usage: ./low_support.sh in.paf ver platform
  in.paf  : bam2paf
  ver     : assembly version prefix
  platform: HiFi or ONT
```
Adjust the file paths below before running this script.

### Prerequisites

Unlike the rest, this script requires several file paths. `asm` and `telo` can be generated with [`init.sh`](init.sh).

* asm    : `ver`.bed; Bed file containing regions as chromosome (scaffold) length.
* telo   : `ver`.telo.bed; Ends of chromosomes to exclude from being flagged as low coverage due to natural sequencing coverage dropouts. Could be an empty file.
* error. : `ver`.error.bed; Regions of known consensus base errors. [Merqury](https://github.com/arangrhie/T2T-Polish/tree/master/merqury) asm only track.
* pattern: pattern/`ver` folder; containing results from [seqrequester](https://github.com/arangrhie/T2T-Polish/tree/master/pattern).

### HiFi mode

Asset runs without excluding any bases at the beginning or end of a read alignment while collecting coverage (`-l 0`). Clipping threshold is set to `10`, flagging regions when > 10 reads or > 10% of the reads have clippings.

### ONT mode

Asset runs with default options. Clipping threshold is set to `15`, flagging regions when > 15 reads or > 15% of the reads have clippings.

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
