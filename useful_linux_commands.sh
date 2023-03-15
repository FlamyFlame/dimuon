# use ${filename%.png} to get rid of the subfix .png

for n in $(ls *png); do newn=${n%.png}; mv $n $newn"_noDPHI_div.png"; done
