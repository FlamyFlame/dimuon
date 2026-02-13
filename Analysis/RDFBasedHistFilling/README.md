# RDataFrame-based Histogram Filling Classes

## Class Structure & Overview
* 

## What User Must Implement
### I/O Variables
* `input_files`: `std::vector<std::string>`	- list of input files where `muon_pair_tree_sign1` & `muon_pair_tree_sign2` live
* `output_file`: `std::string` - output file path
* `infile_var1D_json`: `std::string` - file path for 1D-variable json file

### Pure Virtual function hooks 
* `SetIOPathsHook` - set I/O paths
* `FillHistograms` - fill histograms by building data samples of interest & calling `FillHistogramsSingleDataFrame` on each data sample

### Virtual function extras
* `InitAnalysisSettingsHook`
* `BuildHistBinningMapExtra`
* `BuildFilterToVarListMapExtra`
* `HistPostProcessExtra`
* `WriteOutputExtra`
* `CleanupExtra`

### Notes for Pythia & Powheg
* `meta_tree` entry `nentries_before_cut` is needed for filter efficiencies, where `filter effcy = P[dimuon with 3.7GeV & |eta| < 2.5 cuts | MC hard scattering mode] × P[dimuon with 4GeV & |eta| < 2.4 & additional cuts | dimuon with 3.7GeV & |eta| < 2.5 cuts]` gets multiplied to the hard-scattering cross sections to give the final cross section of dimuon events
	* When `hadd` multiple batches, the entries in `meta_tree` gets concatenated --> loop over & add them if needed.