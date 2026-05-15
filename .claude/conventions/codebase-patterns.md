# Codebase Patterns and Conventions

## ACLiC compilation

- `RDFBasedHistFillingPP.cxx` and `RDFBasedHistFillingPbPb.cxx` both `#include "RDFBasedHistFillingData.cxx"` (no pragma once). Must compile in **separate ROOT sessions** — compiling both in one session causes redefinition errors.
- If `redefinition of 'RDFBasedHistFillingData'` appears, delete stale `*_cxx.so`, `*_cxx.d`, `*_cxx_ACLiC_dict*.pcm` files before retrying.
- Do NOT preload `RDFBasedHistFillingPowheg.cxx` before compiling PP or PbPb.
- When adding a new derived class that includes Data.cxx, give it its own test script that compiles it in isolation.

## DatasetTriggerMap coordination

`DatasetTriggerMap.h` maps (year, data_type) to file-name suffix. Must match `base_trig_suffix` set by `trigger_mode` in `RDFBasedHistFillingData.cxx`. If they disagree, the combined plotter silently skips the year.

## Histogram binning registration

Every named binning string in `var1D_*.json` must be registered in `hist_binning_map` (via `BuildHistBinningMap*()` functions). If a name is missing, `ReadVar1DJson()` crashes during `Initialize()`.

Registration locations:
- `BuildHistBinningMapBaseCommon()` in `RDFBasedHistFillingBaseClass.cxx` — shared binnings
- `BuildHistBinningMapDataCommon()` in `RDFBasedHistFillingData.cxx` — data-specific
- `BuildHistBinningMapPbPbExtra()` in `RDFBasedHistFillingPbPb.cxx` — PbPb-specific

The var1D_dict key is `hist_name` (e.g. `"pair_pt_log"`), not `var_name` (e.g. `"pair_pt"`).

## CRTP mixin architecture

NTuple processing uses CRTP with variadic Extras mixins:
```
DimuonAlgCoreT<PairT, MuonT, Derived, ...Extras>
```
Extras are composed via multiple inheritance (PPExtras, PbPbExtras, PythiaTruthExtras, etc.).
Concrete classes are in `*AnalysisClasses.h`. See `Analysis/README.md` for the full hierarchy.

## Sign conventions

- `muon_pair_tree_sign1` = same-sign (SS), suffix `_ss`
- `muon_pair_tree_sign2` = opposite-sign (OS), suffix `_op`

## TH3DModel (ROOT 6.34)

ROOT 6.34 `TH3DModel` has only 8-arg (all variable-bin) and 11-arg (all uniform) constructors. For mixed axes, generate explicit edge vectors for uniform axes and use the 8-arg form.
