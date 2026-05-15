# RDataFrame Patterns

## Differential scaling

For differential plots (dN/dX, dsigma/dX), always use:
```cpp
h->Scale(NORM_FACTOR, "width");
```
ROOT divides each bin by its own width automatically. Never manually compute bin_width and divide — it breaks for non-uniform binning.

## var1D JSON registration

Histogram variables are defined in `var1D_*.json` files. Each entry specifies `hist_name`, `var_name`, and `binning` (either numeric edges or a named string).

Named binning strings must be registered in `hist_binning_map` before use. See `conventions/codebase-patterns.md` for registration locations.

## Pair eta binning

For pair-level plots (reco efficiency, cross-section, signal acceptance), use `pair_eta_proj_ranges_coarse_incl_gap` from `CommonEffcyConfig.h`.

Do NOT use `q_eta_proj_ranges_*` for pair eta — those are for single-muon trigger efficiency fitting using q*eta only.

## Compilation isolation

PP and PbPb hist-filling classes must be compiled in separate ROOT sessions. See `conventions/codebase-patterns.md` for details.

## IsPbPb() virtual dispatch

`IsPbPb()` is a virtual method (base returns false, PbPb overrides to true) in `RDFBasedHistFillingData.h`. Use this instead of `dynamic_cast` to distinguish PP vs PbPb at runtime in shared code.
