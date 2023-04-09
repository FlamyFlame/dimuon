# use ${filename%.png} to get rid of the subfix .png

for n in $(ls *png); do newn=${n%.png}; mv $n $newn"_noDPHI_div.png"; done
for n in $(ls run_muon_pair_plotting*sh); do newn=$(echo $n | sed 's/run_muon_pair_plotting/run_plotting/g'); mv $n $newn; done
for n in $(ls *_mc_data_compr*); do newn=$(echo $n | sed 's/_mc_data_compr/_powheg_data_compr/g'); mv $n $newn; done
for n in $(ls *_powheg_pythia_data_compr*); do newn=$(echo $n | sed 's/_powheg_pythia_data_compr/_mc_data_compr/g'); mv $n $newn; done
for n in $(ls mc_muon_pair*); do newn=$(echo $n | sed 's/mc/pythia/g'); mv $n $newn; done
for n in $(ls hard_scatt*); do newn=$(echo $n | sed 's/hard_scatt/hard_scatt_near_away/g'); mv $n $newn; done