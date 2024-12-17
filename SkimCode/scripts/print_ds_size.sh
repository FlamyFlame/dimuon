# script to write the list of data23_hi datasets in the physics_HardProbes stream into a file ds_size.txt
# along with the dataset size (in unit of TB/GB/MB)
# should copy into the run directory for data23_hi data before running this script

rm -f ds_size.txt;
touch ds_size.txt;

for file in $(rucio list-dids --short data23_hi.*.physics_HardProbes.merge.AOD.f*); do 
	echo $file >> ds_size.txt; 
	rucio list-files $file | tail -2 | head -1 >> ds_size.txt; 
done