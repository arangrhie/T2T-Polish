#!/bin/sh

:<<'END'
for seq in ONT clr hifi
do
	mkdir -p ${seq}_all
	cd ${seq}_all
	ln -s ../ast.sh
	ln -s /data/Phillippy/t2t-share/team-curation/alignments/20200904/$seq/output.paf $seq.paf
	cd ../

	mkdir -p ${seq}_pri
	cd ${seq}_pri
	ln -s ../ast.sh
	ln -s /data/Phillippy/t2t-share/team-curation/alignments/20200904/$seq/output.primary.paf $seq.primary.paf
	cd ../
done
#END


seq=clr
cd ${seq}_all
ln -s ../ast.sh
./ast.sh ${seq}_all 35
#cd ../${seq}_pri
#./ast.sh ${seq}_pri 35
cd ../
END

#:<<'END'
seq=hifi
cd ${seq}_all
#ln -s ../ast.sh
#./ast.sh ${seq}_all 31
cd ../${seq}_pri
ln -s ../ast.sh
# Make the mean cov to something high so we get 'reliable' regions for everywhere >10 reads
# Setting this to 31 will generate the normal reliable tracks, where regions with mean_cov*2.5 gets flagged with the <10 read region
./ast.sh ${seq}_pri 100000000
cd ../
#END

:<<'END'
seq=ONT
cd ${seq}_all
#./ast.sh ${seq}_all 127
cd ../${seq}_pri
#./ast.sh ${seq}_pri
./ast.sh ${seq}_pri 127
cd ../
END
