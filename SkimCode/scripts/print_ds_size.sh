rm -f ds_size.txt;
touch ds_size.txt;

for file in $(rucio list-dids --short data23_hi.*.physics_HardProbes.merge.AOD.f*); do 
	echo $file >> ds_size.txt; 
	rucio list-files $file | tail -2 | head -1 >> ds_size.txt; 
done