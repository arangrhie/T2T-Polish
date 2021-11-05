# Microsatellites

Collect 2-mer microsatellites (simple tandem repeats) in every 128 bp window.

## Dependencies
* [seqrequester](https://github.com/marbl/seqrequester)
* java

```
Usage: ./microsatellites.sh in.fasta
```

Output:
* 4 wig files: GA, TC, GC, and AT
* 4 bed files: Regions over 80%

`seqrequester microsatellite` collects pattern of sequences (`ga` `gc` or `at`) in bed files.
`bedToFixedWig.jar` converts the `.bed` file to `.wig` file.
Additionally, regions with > 80% are collected with `gt80.bed` postfix.


The `.gt80.bed` files are used for generating the `issues.bed` files to find
[low coverage regions associated with sequencing biases](https://github.com/arangrhie/T2T-Polish/tree/master/coverage#prerequisites-1).
