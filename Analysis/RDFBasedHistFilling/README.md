# RDataFrame-based Histogram Filling Classes

## Class Structure & Workflow
### Initialize
* Initialize I/O & analysis-specific settings - Must be specified by child analysis classes
* `BuildFilterToVarListMap`: Map data-frame/filter describer strings to 1/2/3D variables to be plotted

* `BuildHistBinningMap`: Map histogram-binning describer strings to actual binnings
* `ReadVar1DJson`: Read json configuratin file (`infile_var1D_json`) for dictionary of 1D variables

### ProcessData
* `CreateBaseRDFs` & `FillHistograms`: Create RDataFrames & fill histograms
	* Call function `FillHistogramsSingleDataFrame` for each data frame (data sample)
	* `FillHistogramsSingleDataFrame`: 
		1. retrieves list of 1/2/3D variables from either filter or (filter, weight) pairs
		2. Calls `Histo1/2/3D` & pushes result into `hist1/2/3d_rresultptr_map`
		* Include options to plot weighted histograms & to NOT write output 1/2/3D histograms
* `HistPostProcess`: histogram post processing
	* Runs RDF analysis & cast `ROOT::RDF::RResultPtr< ::TH1/2/3D>>` in `hist1/2/3d_rresultptr_map` to `TH1/2/3D*` & store in `hist1/2/3D_map`
	* Perform user-specified post processing: addition / division / projection / summary-statistics calculation etc.

### Finalize
* `WriteOutput`: Write outputs: all histograms in `hist1/2/3D_map` except for those listed in `hists_to_not_write` & user-specified extras
* `Cleanup`: Close output file & delete raw pointers

## What User Must Implement
### I/O Variables
* `input_files`: `std::vector<std::string>`	- list of input files where `muon_pair_tree_sign1` & `muon_pair_tree_sign2` live
* `output_file`: `std::string` - output file path
* `infile_var1D_json`: `std::string` - file path for 1D-variable json file

### Key Variables for Histogram Plotting - MUST BE IMPLEMENTED BY USERS
* `df_filter_to_var1/2/3D_list_map` & `df_filter_and_weight_to_var1/2/3D_list_map`: Map of data-frame/filter describer strings to 1/2/3D variables to be plotted for this data frame
	* Implement in `BuildFilterToVarListMap`
* `hist_binning_map`: Map histogram-binning describer strings to actual binnings
	* Implement in ``
* `weight_specifier_to_column_map`: map of weight-specifier string to column name to be used for weight
	* Use if weights involved
* `df_map`: Map of data-frame describer strings to data frames - index for data frames, useful for searching
	* Implement & used in `CreateBaseRDFs` & `FillHistograms`

### Pure Virtual function hooks - MUST BE IMPLEMENTED BY USERS
* `SetIOPathsHook` - set I/O paths
* `FillHistograms` - fill histograms by building data samples of interest & calling `FillHistogramsSingleDataFrame` on each data sample

### Virtual function extras - should be implemented/checked by users
* `InitAnalysisSettingsHook` 		- user/analysis-specific settings
* `BuildHistBinningMapExtra` 		- user-specific histogram binning specifications
* `BuildFilterToVarListMapExtra` 	- user-specific filter to variable list maps
* `HistPostProcessExtra` 			- user-specific histogram post processing
* `WriteOutputExtra` 				- user-specific output writings
* `CleanupExtra`					- user-specific cleanups

### Notes for Pythia & Powheg
* `meta_tree` entry `nentries_before_cut` is needed for filter efficiencies, where `filter effcy = P[dimuon with 3.7GeV & |eta| < 2.5 cuts | MC hard scattering mode] × P[dimuon with 4GeV & |eta| < 2.4 & additional cuts | dimuon with 3.7GeV & |eta| < 2.5 cuts]` gets multiplied to the hard-scattering cross sections to give the final cross section of dimuon events
	* When `hadd` multiple batches, the entries in `meta_tree` gets concatenated --> loop over & add them to use for weight.

