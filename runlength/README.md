# Homopolymer and 2-mer microsatellite run length matrix

The run-length matrix using `runLengthMatrix` executable is available with [Margin](https://github.com/UCSC-nanopore-cgl/margin).

```bash
./runLengthMatrix \
<BAM> \
<REF> \
base_params.json \ 
-t <THREADS> \
-l <highest_runlength_allowed> \
-o OUTPUT_DIR/OUTPUT_PREFIX

base_params.json = https://github.com/UCSC-nanopore-cgl/margin/blob/master/params/base_params.json
```

