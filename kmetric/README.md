# K*metric
K* related scripts.

## Useful oneliners
Create .Wig file from K*
```
awk 'BEGIN{print "track autoScale=on"}{if($1!=chr){chr=$1; counter=1; print "variableStep chrom="chr" span=1"};if($6=="NA"){printf counter"\t"0}else{print counter"\t"$6}counter+=1}' t2t-chm13.20200727.asm.combined
```
