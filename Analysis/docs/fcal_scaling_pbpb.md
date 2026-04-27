# FCal ET Reweighting for PbPb 2024/2025 Centrality

## Motivation

PbPb 2024 and 2025 do not yet have Glauber-calibrated FCal centrality thresholds.
The PbPb 2023 thresholds (`FCal_ET_Bins_PbPb2023`) are used as reference.
To apply them consistently, the FCal ET distribution of 2024/2025 is reweighted to match 2023
using a per-event, FCal-ET-dependent correction weight.

## Procedure

### Step 1 — Derive per-year FCal weights (`plotting_codes/fcal_scaling/derive_fcal_scaling.cxx`)

1. For each year (23, 24, 25): apply the full 5-cut event selection to the raw skim using the
   per-year cut file (`event_sel_cuts_pbpb_20YY.root`).
2. Fill FCal ET histogram: `FCal_Et = (FCal_Et_P + FCal_Et_N) * 1e-6` [TeV].
3. Restrict to centrality 0-80% region: `FCal_Et > kFCal_80_ref = 0.063208 TeV`
   (= `FCal_ET_Bins_PbPb2023[79]`, the boundary between 79% and 80% centrality in 2023).
4. Normalise each year's histogram to area = 1 in this region.
5. Per-bin weight: `w(FCal_bin) = h23_norm(FCal_bin) / hyr_norm(FCal_bin)`.
   Bins where `h_yr = 0` are skipped (TGraph::Eval interpolates across them).
   Weights are capped at 5 to suppress statistical noise at the high-FCal tail.
6. Save `TGraph("g_fcal_weight")` to `pbpb_20YY/fcal_weight_pbpb_20YY.root` for years 24, 25.
   Year 23 is the reference; no file needed (weight = 1 everywhere).
7. Output plots to `plots/fcal_scaling/`:
   - `fcal_before_scaling.png`: 2-panel raw comparison (2024/2025 vs 2023).
   - `fcal_after_scaling.png`: 2-column ratio plot — upper panels: weighted FCal distributions
     (normalised over FCal > 80% boundary) overlaid on 2023; lower panels: weight(FCal) per year.

### Step 2 — Load weight in NTuple processing (`NTupleProcessingCode/PbPbExtras.{h,c}`)

At initialisation (`InitEventSel`), for years 24 and 25:
load `TGraph("g_fcal_weight")` from the weight file; store as member `g_fcal_weight_`
(nullptr for year 23).

At fill time (`FillMuonPairExtra`): raw FCal ET is stored in the pair unchanged.
The per-event weight is evaluated and stored in `PairPbPbExtras::fcal_corr_weight`:

```
FCal_Et > kFCal_80_ref (centrality 0-80%): fcal_corr_weight = g_fcal_weight_->Eval(FCal_Et)
FCal_Et ≤ kFCal_80_ref (centrality >80%):  fcal_corr_weight = 0
Year 23:                                    fcal_corr_weight = 1 (default, no file loaded)
```

### Step 3 — Centrality recalculation (`MuonObjectsParamsAndHelpers/MuonPairPbPb.h`)

For years 24 and 25, `PairValueCalcHook` calls `UpdateCentrality()`, which applies the
2023 FCal thresholds (`GetCentralityPbPb2023`) to the raw FCal ET.

`UpdateCentrality()` is called from `MuonPairBase::Update()` → `PairValueCalc()` →
`PairValueCalcHook()`, which runs **after** `FillMuonPairExtra()` in the event loop
(see `DimuonDataAlgCoreT.c`). So `fcal_corr_weight` is always set before centrality
is recalculated.

### Step 4 — Application in histogram filling (`RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx`)

`fcal_corr_weight` is stored as a branch in the output NTuple and must be applied
when filling physics histograms. It is folded into `weight_for_RAA` in
`FillHistogramsCrossx`:

```cpp
df.Define("weight_for_RAA",
    [this](int avg_centrality, double weight, float fcal_corr_weight) {
        return this->CalculateWeightForRAA(avg_centrality, weight) * fcal_corr_weight;
    },
    {"avg_centrality", "weight", "fcal_corr_weight"});
```

`CalculateWeightForRAA` contributes the per-centrality-bin T_AA / lumi factor;
`fcal_corr_weight` corrects the FCal distribution shape.

For year 23, `fcal_corr_weight = 1.0` (no weight file loaded) → no change.
For events with FCal ≤ 80% boundary, `fcal_corr_weight = 0` → those events are
automatically excluded from all weighted histograms.

## Per-event flow summary

| Year | Event selection | FCal stored | Centrality | `fcal_corr_weight` |
|------|----------------|-------------|------------|---------------------|
| 23   | raw FCal, per-23 cuts | raw | from `centrality` branch | 1.0 |
| 24   | raw FCal, per-24 cuts | raw | recalc via `GetCentralityPbPb2023` | `g_fcal_weight_->Eval(FCal)` for FCal > boundary; 0 otherwise |
| 25   | raw FCal, per-25 cuts | raw | recalc via `GetCentralityPbPb2023` | same |

## Key constants and files

| Quantity | Value / Location |
|----------|-----------------|
| 80% centrality boundary | `kFCal_80_ref = 0.063208 TeV` = `FCal_ET_Bins_PbPb2023[79]` |
| Weight TGraph key | `"g_fcal_weight"` in `fcal_weight_pbpb_20YY.root` |
| Path function | `PbPbFCalWeightPath(run_year)` in `PbPbEventSelConfig.h` |
| Weight member | `TGraph* g_fcal_weight_` in `PbPbExtras` |
| Per-event weight field | `float fcal_corr_weight` in `PairPbPbExtras` |

## Impact on final results

**Differential cross-section**: within each centrality bin, events are weighted by
`fcal_corr_weight`. This corrects the FCal composition inside the bin to match 2023,
so the geometric interpretation (T_AA, ⟨N_coll⟩) from 2023 Glauber is valid.

**RAA**: `RAA = (weighted PbPb yield / T_AA) / (pp reference)`. T_AA from 2023 Glauber
is self-consistent with the centrality selection (same thresholds). `fcal_corr_weight`
corrects the FCal distribution shape; it does **not** correct overall normalization
(that is handled by luminosity / N_inel).

**Trigger efficiency histograms**: filled without `fcal_corr_weight` (unweighted counts),
which is correct — efficiency ratios are independent of the FCal reweighting.

**Statistical cost**: the weight redistributes effective statistics between FCal bins.
The effect is small when 24/25 FCal scales are close to 2023.

## Notes

- The event selection always uses raw FCal ET — per-year cut curves were derived on the raw scale.
- `fcal_corr_weight` is applied in `FillHistogramsCrossx` (both the global and
  per-centrality-bin nodes). Generic trigger-efficiency histograms are unweighted.
- Centrality bin assignment uses raw FCal with 2023 thresholds, accepting the small systematic
  from any residual FCal scale difference between years (~few percent).
- If proper Glauber calibrations become available for 2024/2025: update `FCal_ET_Bins_PbPb2024/2025`,
  remove the `UpdateCentrality()` calls for those years, and retire the weight files.
