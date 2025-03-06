# Directory for dimuon-related code

## Analysis 
**contains all pieces of code for AFTER-ntuple analysis, including**
	- ROOT file of event-based TTrees --> ROOT file of muon-pair-based TTrees
	- histogram making: muon-pair-based TTrees --> ROOT file containing all histograms
	- plotting: ROOT file of histograms --> PDF/PNG files with campus division, plotting style, etc., executed


## SkimCode 
**athena-based packages for skimming data in the format of AOD files and produce ntuples, including specifically**
	- source directory --> source code of the packages themselves
	- CMakeFiles
	- scripts for python configuration files + grid-submitting bash scripts
	- XMLs: including GRLs for different data periods
	- datasetnames



## Paths for data files, muon pairs, hist root file, and final output plots (png files)
	- data
		- directory: `/usatlas/u/yuhanguo/usatlasdata/dimuon_data`
			- subdirectory for data files, muon-pair trees & hist root files: pp_run2, pp_2024, pbpb_run2, pbpb_2023, pbpb_2024
			- subdirectory for data & data-MC-comparison plots: plots/
	- pythia
		- directory: `/usatlas/u/yuhanguo/usatlasdata/pythia`
		- pythia sample files, muon-pair trees & hist root files are at this level
		- subdirectory for pythia plots: plots/
	- powheg
		- `/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample` --> bb_full_sample, cc_full_sample
			- powheg sample files, muon-pair trees & hist root files are at this level in the bb_full_sample, cc_full_sample subdirectories
		- powheg plots are at `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/powheg/`
