bedtools makewindows -w 1000000 -g t2t-chm13.20200904.fasta.fai > 1M.windows.bed
bedtools coverage -a 1M.windows.bed -b sv_small_no_rDNA.vcf > 1M.windows.dv.cov.bed
bedtools coverage -a 1M.windows.bed -b racon.round3.vcf > 1M.windows.rc.cov.bed
bedtools coverage -a 1M.windows.bed -b racon.round3.merfin.vcf > 1M.windows.rm.cov.bed
