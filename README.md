# T2T-Polish

This repository contains up-to-date evaluation and polishing workflows to adapt on general genome assembly projects, with most of the ideas developed and described in [this paper](https://doi.org/10.1101/2021.07.02.450803).

For exact command lines and workflows used to generate the T2T-CHM13v1.0 and T2T-CHM13v1.1 assemblies, please refer to the [Methods](https://github.com/marbl/CHM13-issues#methods) section in the [CHM13-Issues](https://github.com/marbl/CHM13-issues) repo. Note that some of the tools have been updated since then, and are tracked on this repo.

## Contents
* [QV estimate with hybrid k-mer db](merqury)
* [Homopolymer and 2-mer microsatellite run length comparison](runlength)
* [K\* metric](kmetric)
* Repeat-aware [Winnowmap2 alignments](winnowmap) and [Marker assisted alignment filtering](marker_assisted)
* [Coverage analysis](coverage)
* [Polish SV and SNV like errors through a case study on Chr. 20](doc/T2T_polishing_case_study.md)
* [Automated polishing with Racon + Merfin](automated_polishing)

## Related external links

### Variant call, refinements and formatting (Also see [Error Detection](https://github.com/marbl/CHM13-issues/blob/main/error_detection.md))
* [Call SVs and refine sequences](https://github.com/malonge/CallSV)
* [Merge vcfs from Illumina/HiFi and ONT](https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/polishing_merge_script/vcf_merge_t2t.py)
* [Generate telomere edits](https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/telomere_variants/generate_telomere_edits.py)

### Repeat-aware alignments
* [Winnowmap](https://github.com/marbl/Winnowmap)
* [TandemMapper2 (re-named as VerityMap)](https://github.com/ablab/VerityMap)

### Automated polishing
* [Racon](https://github.com/isovic/racon/tree/liftover): Liftover branch for outputting edits in `.vcf`
* [Merfin](https://github.com/arangrhie/merfin): Latest stable code-base

### Base level QV estimation
* [Meryl](https://github.com/marbl/meryl)
* [Merqury](https://github.com/marbl/merqury)

## Citation
Please cite if any of the codes shared in this repo was used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. bioRxiv (2021) doi: https://doi.org/10.1101/2021.07.02.450803
