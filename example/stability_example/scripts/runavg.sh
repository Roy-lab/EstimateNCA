
       
lst=lists/all_files.txt
tfa=tfa/tfa0.010.txt
echo python2.7 scripts/takeAvg.py --tfalist ${lst} --tfaout ${tfa}
python2.7 scripts/takeAvg.py --tfalist ${lst} --tfaout ${tfa}
for ((k=0;k<10;k++)); 
do 
	let i=k*10;
	let j=i+9;
	lst=lists/files_${i}_to_${j}.txt
	tfa=tfa/tfa0.010_${k}.txt; 
	echo python2.7 scripts/takeAvg.py --tfalist ${lst} --tfaout ${tfa} ; 
	python2.7 scripts/takeAvg.py --tfalist ${lst} --tfaout ${tfa} ; 
done; 
