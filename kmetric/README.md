# K*metric
K* related scripts.

## Useful one liners
Create .Wig file from K*
```
awk 'BEGIN{print "track autoScale=on"}{if($1!=chr){chr=$1; counter=1; print "variableStep chrom="chr" span=1"};if($6=="NA"){printf counter"\t"0"\n"}else{printf counter"\t"$6"\n"}counter+=1}' t2t-chm13.20200727.asm.combined > t2t-chm13.20200727.Wig
```
From merfin dump output (0-based) including missing (note: Wig is 1-based):
```
awk 'BEGIN{print "track autoScale=on"}{if($1!=chr){chr=$1; print "variableStep chrom="chr" span=1"};if($3!=0){printf $2+1"\t"$5"\n"}}' t2t-chm13.20200904 > t2t-chm13.20200904.illumina.Wig```
```
Convert to bigWig
```
wigToBigWig t2t-chm13.20200727.Wig chrs.size t2t-chm13.20200727.kstar.bw
```