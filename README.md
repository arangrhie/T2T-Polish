# T2T-Polish
This repository is a shared code base from the polishing team, created over the T2T Workshop (Summer, 2020) and described in [this paper](https://doi.org/10.1101/2021.07.02.450803).

* [winnowmap](https://github.com/arangrhie/T2T-Polish/tree/master/winnowmap): Launcher script to run Winnowmap in parallel and filter for primary alignments
* [coverage](https://github.com/arangrhie/T2T-Polish/tree/master/coverage): Coverage based analysis for detecting low-coverage regions or regions with excessive clippings
* [marker_assisted](https://github.com/arangrhie/T2T-Polish/tree/master/marker_assisted): Filter mapping file (bam/cram) to generate alignments that have support from markers
* [merfin](https://github.com/arangrhie/T2T-Polish/tree/master/merfin): Scripts we used to run Merfin for polishing T2T-CHM13v0.9.
 Latest code is available [here](https://github.com/arangrhie/merfin).
* [variant_call](https://github.com/arangrhie/T2T-Polish/tree/master/variant_call): Scripts to call variants

## Related external links
* [homopolymer and microsatellite run length](https://github.com/UCSC-nanopore-cgl/margin)
* [script to merge vcfs from Illumina/HiFi and ONT](https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/telomere_variants/generate_telomere_edits.py)

### Repeat aware alignment
* [Winnowmap](https://github.com/marbl/Winnowmap)
* [TandemMapper](https://github.com/ablab/TandemTools)

### Automated polishing
* [Racon](https://github.com/isovic/racon/tree/liftover): Liftover branch for outputting .vcf
* [Merfin](https://github.com/arangrhie/merfin): Latest stable code-base

### Base level QV estimation
* [Merqury](https://github.com/marbl/merqury)

## Citation
Please cite if any of the codes shared in this repo was used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. bioRxiv (2021) doi: https://doi.org/10.1101/2021.07.02.450803
