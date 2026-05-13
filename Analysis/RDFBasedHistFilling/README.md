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

## JSON & Source Code Coordination

### How `var1D_json` binning works
`ReadVar1DJson()` resolves `"binning"` fields in one of two ways:

1. **Inline object** `{"nbins": N, "min": lo, "max": hi}` — uniform bins, no source change needed.
2. **Named string** `"pT_bins_120"` — looks up `hist_binning_map["pT_bins_120"]`.

**Critical**: every named binning used in any JSON must be registered in `hist_binning_map`.  The
base class `BuildHistBinningMapBaseCommon()` registers `pT_bins_40`, `pT_bins_80`, and
`pT_bins_120` (log-spaced, 15 bins 8–120 GeV).  Dataset-specific binnings (e.g.
`pT_bins_single_muon`, `eta_bins_trig_effcy`, centrality-dependent bins) are added via
`BuildHistBinningMapDataCommon()` / `BuildHistBinningMapPbPbExtra()`.  If you add a new named
binning to a JSON you *must* also register it in the appropriate `Build*` function; otherwise
`ReadVar1DJson()` throws and the run crashes before any histograms are filled.

### `var1D_dict` key is `hist_name`, not `var_name`
`var1D_dict` is keyed by the JSON `"hist_name"` field.  The variable looked up in
`df_filter_(and_weight_)to_var1D_list_map` must match `hist_name`, not `var_name`.  Example: the
entry with `"hist_name": "pair_pt_log"` and `"var_name": "pair_pt"` stores its bins under the
key `"pair_pt_log"`.

### Hardcoded 2D/3D histograms
Crossx 2D/3D histograms (e.g. `h2d_crossx_pair_pt_pair_eta_*`) are built directly in
`FillHistogramsCrossx()` using `pms.pT_bins_120.data()` — they bypass the JSON mechanism
entirely.

### `TH3DModel` constructor limitation (ROOT ≥ 6.34)
ROOT 6.34 `TH3DModel` only has an **11-arg all-uniform** and an **8-arg all-variable-bin**
constructor; there is no mixed form.  For axes with log/custom bins on one dimension and uniform
on the others, generate explicit edge vectors and use the 8-arg constructor:
```cpp
auto make_edges = [](int n, double lo, double hi){ ... };
const auto eta_edges = make_edges(44, -2.4, 2.4);
ROOT::RDF::TH3DModel(name, title, npt, ptbins, 44, eta_edges.data(), 50, z_edges.data())
```

## ACLiC Compilation Order

`RDFBasedHistFillingPP.cxx` and `RDFBasedHistFillingPbPb.cxx` both `#include
"RDFBasedHistFillingData.cxx"`.  Because `Data.cxx` lacks `#pragma once`, compiling both in the
**same ROOT session** via `.L PP.cxx+` followed by `.L PbPb.cxx+` (or vice versa) causes ACLiC's
dictionary umbrella to include `Data.cxx` twice → redefinition errors.

Rules:
- **Always compile PP and PbPb in separate ROOT sessions** (i.e. separate `root` invocations).
- `IsPbPb()` is a virtual method (base returns `false`, PbPb overrides to `true`) that replaces
  `dynamic_cast<RDFBasedHistFillingPbPb*>` in `FillHistograms()`, eliminating the inter-class
  typeinfo dependency at link time.
- Do not preload unrelated `.cxx` files (e.g. `RDFBasedHistFillingPowheg.cxx`) before
  compiling PP or PbPb; it is unnecessary and can introduce stale-library conflicts.
- If you see `redefinition of 'RDFBasedHistFillingData'` during compilation, delete stale
  `*_cxx.so`, `*_cxx.d`, and `*_cxx_ACLiC_dict*.pcm` files before retrying.

## Running crossx filling & plotting

```bash
# Fill histograms (PbPb before PP, since run_all_crossx.sh now enforces this order)
cd RDFBasedHistFilling
bash run_crossx_hist_filling_pbpb23.sh
bash run_crossx_hist_filling_pbpb24.sh
bash run_crossx_hist_filling_pbpb25.sh
bash test_crossx_pp24.sh

# Plot
cd ../Analysis
root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx'
root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx(23)'
root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx(24)'

# Or run everything at once:
bash plotting_codes/single_b_analysis/run_all_crossx.sh
```

### Notes for Pythia & Powheg
* `meta_tree` entry `nentries_before_cut` is needed for filter efficiencies, where `filter effcy = P[dimuon with 3.7GeV & |eta| < 2.5 cuts | MC hard scattering mode] × P[dimuon with 4GeV & |eta| < 2.4 & additional cuts | dimuon with 3.7GeV & |eta| < 2.5 cuts]` gets multiplied to the hard-scattering cross sections to give the final cross section of dimuon events
	* When `hadd` multiple batches, the entries in `meta_tree` gets concatenated --> loop over & add them to use for weight.

