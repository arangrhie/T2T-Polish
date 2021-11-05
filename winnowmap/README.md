# Winnowmap Alignment
Align long-reads with [Winnowmap](https://github.com/marbl/Winnowmap/releases).

From the reads, top 0.02% of the repetitive 15-mers were collected using Meryl and downweighted in the process of minimizer sampling.
Use platform specific mapping parameters;
`map-pb` for HiFi, `map-pb-clr` for CLR, and `map-ont` for ONT reads, respectively.
Reads are sorted and indexed using samtools, and filtered for primary reads with `samtools view -F0x100`.

## 1. init.sh
```shell
# Collect repetitive 15-mers
meryl count k=15 $ref output merylDB
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
```

## 2. map.sh
This step can be run in parallel per input read file. Use [`split.sh`](split.sh) to split a file by every X num. lines (e.g. 100 M).

```shell
# Map and sort reads
winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t$cpus $ref $reads > $out.sam
samtools sort -@$cpus -m2G -T $out.tmp -O bam -o $out.sort.bam $out.sam
```

## 3. merge.sh
```shell
# Merge and index
samtools merge -O bam -@$SLURM_CPUS_PER_TASK $prefix.bam *.sort.bam
samtools index $prefix.bam
```

## 4. filter.sh
```shell
# Filter to get primary read alignments
samtools view -F0x104 -@$SLURM_CPUS_PER_TASK -hb $in > $out
samtools index $out
```

## 5. sam2paf.sh
Finally, convert the desired bam file to paf using [paftools from minimap2](https://github.com/lh3/minimap2/tree/master/misc)
and processe to coverage analysis. This script can be found under [coverage](https://github.com/arangrhie/T2T-Polish/tree/master/coverage/sam2paf.sh).
```shell
./sam2paf.sh $prefix.pri.bam $prefix.pri.paf
```

Optionally, spurious alignment blocks shorter than 1 kb or lower than 85% identity can be removed from the ONT primary alignments.

```shell
# Filter ONT for alignment blocks < 1kb and identity < 85%
awk '$11>1000 && $10/$11>0.85' ont.pri.paf > ont.pri.len1k_idy85.paf
```

## On Slurm
The entire job can be submitted using `_submit.sh`.
```shell
Usage: ./_submit.sh ref.fasta prefix map [wm_opt]
  Required: input.fofn
```

Collect all input reads to input.fofn, provide `ref.fasta`, `prefix`, and `map` option (`map-pb` / `map-pb-clr` / `map-ont`).
Additional winnowmap options could be provided in `wm_opt`, with quotes.
