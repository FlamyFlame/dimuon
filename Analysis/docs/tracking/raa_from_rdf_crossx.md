# R_AA from RDF crossx (modernize RAA_plotting; consolidate SingleBAnalysis → RDF)

**Mode:** Implementation. **Created:** 2026-06-16. **Roadmap:** task_06.
**Reviewer rules:** code → `/review-analysis-code`; plots → `/review-plot`.
**Depends on:** reco-eff placeholder applied to nominal crossx (DONE,
`reco_eff_placeholder_run2.md`). **Sibling:** `analysis_overview.md` (§2,§4),
`placeholder.md` (items 1–3), `roadmap` Q2.3 (T_AA), `task_06_raa.md`.

---

## Objective
Produce the **reco-corrected R_AA** from the current RDF crossx outputs.
`RAA_plotting.cxx` is stale legacy wired to the retired `SingleBAnalysis`
analyzers; modernize it to read the RDF `histograms_real_pairs_*` outputs, and
add the missing **same-sign (SS)** signal-region histograms to the RDF crossx so
the methodology's OS−SS combinatorial subtraction works. Long-term goal (user):
RDF becomes the ONLY hist-filling code; `SingleBAnalysis/` is reference only.

## Physics Procedure (AUTHORITATIVE — from analysis_overview.md §2,§4)

### 1. Motivation
R_AA quantifies b-quark energy loss in QGP: PbPb yield per binary collision
relative to the pp cross-section. Deviation from 1 ⇒ nuclear modification.

### 2. Top-level equation
```
R_AA(X) = (1/⟨T_AA⟩) · (dN_PbPb/dX) / (dσ_pp/dX),   X ∈ {pair p_T, pair η}, per centrality
```
- **dN_PbPb/dX**: per-event-normalized PbPb yield = (T_AA·L)-weighted OS−SS
  yield in the signal region (the crossx histograms; `weight_for_RAA` carries
  1/(σ·T_AA·L) normalization, see `CalculateWeightForRAA`).
- **dσ_pp/dX**: pp differential cross-section = (1/L_pp)-weighted OS−SS yield.
- **⟨T_AA⟩**: Glauber nuclear overlap per centrality. **2023 placeholder** for
  all years (roadmap Q2.3); folded into `weight_for_RAA`/crossx_factor already.
- Both pp and PbPb yields are efficiency-corrected (trigger × reco placeholder).

### 3. Step-by-step
a. **Yield & combinatorial background (§4a):** count OS pairs in the signal
   region (minv∈[1.08,2.9], pair p_T>8, |q·η|<2.2 per muon, dR>0.05). Subtract
   the uncorrelated combinatorial component estimated from **SS** pairs passing
   the *same* cuts: `dN/dX = OS − SS`, bin by bin, BEFORE forming R_AA.
b. **Efficiency:** OS and SS both efficiency-corrected with the same per-muon
   weights (trigger × reco placeholder), so the subtraction is of corrected yields.
c. **Normalization:** PbPb via crossx_factor = 1/(σ_PbPb·f_ctr·⟨T_AA⟩·L_PbPb)
   (per centrality, in `weight_for_RAA`); pp via 1/L_pp (`crossx_weight`).
d. **R_AA:** divide normalized PbPb by normalized pp per (X, centrality), then
   multiply by 1/⟨T_AA⟩ where not already in the weight (it is, via crossx_factor
   — confirm no double T_AA).

### 4. Negative constraints
- MUST subtract SS (combinatorial) from OS in BOTH pp and PbPb before R_AA;
  R_AA of raw-OS would be biased by the uncorrelated background.
- MUST NOT double-count ⟨T_AA⟩ (it is inside the PbPb crossx_factor /
  weight_for_RAA; R_AA must not multiply by 1/T_AA again).
- SS and OS must use IDENTICAL binning, cuts, and efficiency weights.
- pp and PbPb yields must use the SAME pair-pT and pair-η binning for the ratio.
- ⟨T_AA⟩ is a 2023 placeholder — R_AA is preliminary; disclose in the note.

## Context
- Legacy (to consolidate/retire): `SingleBAnalysis/simplified_single_b_analysis_{pp,PbPb}.cxx`
  — read `muon_pairs_*`, write `*_single_b_ana_hists*.root` with
  `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts` (pp) and
  `h3d_op_crossx_w_signal_cuts_vs_centr_vs_pair_eta_vs_pair_pt` (PbPb) — **OS only**,
  binning `pT_bins_80`, `eta_bins`.
- Current RDF crossx (`RDFBasedHistFillingPP/PbPb::FillHistogramsCrossx`) already
  produces those SAME histogram names (now reco+trig corrected), binning
  `pT_bins_120` + 44 η bins + nCtrBins centrality. **OS only** — no SS.
- `RAA_plotting.cxx`: reads `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts`
  (pp) and `h3d_{op,ss}_crossx_..._vs_centr...` (PbPb), does `OS.Add(SS,-1)`;
  hardcoded Mac `base_dir`, legacy filenames, rebins assume `pT_bins_80`.

## Design Decisions
- **SS in RDF, not a new analyzer:** add SS signal-region histograms in
  `FillHistogramsCrossx` mirroring the OS ones (same cuts/weights, `df_ss`),
  rather than reviving SingleBAnalysis. Advances the "RDF-only" goal.
- **Binning:** R_AA uses the RDF binning (`pT_bins_120`, nCtrBins). Rewrite
  RAA_plotting's hardcoded `pT_bins_80`-based rebins to the RDF axes (read bin
  counts from the histograms; group by value ranges, not fixed indices).
- **2025 lumi** already correct (2.59933 nb⁻¹, verified vs CSV); **T_AA = 2023
  placeholder** kept (per user).

## How RAA_plotting.cxx will be modified (per user request)
1. `base_dir`: Mac path → cluster `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/`.
2. Input files: `pbpb_20YY/histograms_real_pairs_pbpb_20YY_single_mu4_no_trg_plots_nominal.root`
   and `pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root` (RDF outputs),
   replacing the `*_single_b_ana_hists*` legacy names. PbPb **combined over years**
   (sum the 3 year files) per analysis convention.
3. Histograms: OS `h3d_op_crossx_..._vs_centr...` + new SS
   `h3d_ss_crossx_..._vs_centr...` (PbPb); OS `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts`
   + new SS `h2d_ss_crossx_pair_pt_pair_eta_binned_w_signal_cuts` (pp). OS−SS.
4. Binning/rebins: derive bin groupings from the actual axes (`pT_bins_120`,
   nCtrBins), not hardcoded `pT_bins_80` indices.
5. Output dir: `dimuon_data/plots/single_b_analysis/RAA/` (cluster). PNG only.
6. SS subtraction guarded: if SS histogram absent, fall back to OS-only with a
   printed WARNING (so the code is runnable during staged rollout).
7. Keep `/review-analysis-code` then `/review-plot`.

## Implementation Plan
1. Tracking doc + Physics Procedure (this file). DONE.
2. Add SS signal-region crossx histograms to RDF (PP + PbPb), same weights as
   OS. ACLiC. → `/review-analysis-code` (quote §3a,§4). 
3. Rerun crossx RDF (pp + 3 PbPb, skip ntuple) → SS histos populated.
4. Modernize `RAA_plotting.cxx` per §"How…above". → `/review-analysis-code`.
5. Run R_AA (modes 1=vs pT, 2=vs η, 3=vs ctr); → `/review-plot`.
6. Update docs (placeholder.md R_AA note, roadmap task_06, status), memory if
   any durable fact.

## Progress Log
- 2026-06-16 — Steps 2–5 DONE.
  - **Step 2** SS histos added to RDF (PbPb 3D `h3d_ss_..._vs_centr`, pp 2D
    `h2d_ss_...`), same cuts/weights as OS. `/review-analysis-code` PASS. ACLiC clean.
  - **Step 3** crossx RDF reran (pp + 3 PbPb, rc=0). SS/OS = 8.8% (pp), 28% (PbPb23) —
    combinatorial larger in PbPb, as expected.
  - **Step 4** `RAA_plotting.cxx` modernized: cluster `base_dir`; new case 6 (Run 3
    PbPb 23+24+25 combined + pp24 2mu4, RDF inputs); OS−SS for pp (2D) and PbPb (3D);
    pp+PbPb both combinatorial-subtracted; member pointers init nullptr (was a
    segfault — uninitialized `h3d_*` used in the combine loop); `pair_pt_rebins`
    updated to 15 bins (RDF pT axis 8–120, was pT_bins_80/12); explicit index-wise
    R_AA ratio (modes 1,2) to avoid the benign TH1::Divide "different bin limits"
    warning (PbPb eta = explicit edge array vs pp eta = uniform; identical edges);
    output `dimuon_data/plots/single_b_analysis/RAA/`; SS-missing guard (OS-only fallback).
  - **Step 5** R_AA ran clean (0 warnings), 3 PNGs (vs pair_pt, pair_eta, ctr).
  - **⚠ NORMALIZATION CAVEAT (surfaced to user):** the combined-year R_AA naive-sums
    the per-year `1/L_year`-weighted 3D histos (matching the EXISTING crossx combined
    plotter `SingleBCrossxPlotterPbPbCombined::GetHistObject`, which also `Add`s).
    This inflates the absolute combined normalization by ~Σ(years) (≈3×): correct
    combined is `ΣN/ΣL`, not `Σ(N_year/L_year)`. Observed R_AA ~1.5–2.5 at mid-pT
    (÷3 → ~0.5–0.8, physical for b-quarks). This is a PRE-EXISTING convention shared
    with the crossx combined plots, NOT introduced by R_AA. Kept R_AA consistent with
    crossx; flagged for a consistent fix (combine raw yields / total lumi) across
    BOTH crossx and R_AA. Absolute scale also preliminary due to 2023 T_AA placeholder.
- 2026-06-16 — Step 1: doc created. Confirmed methodology requires OS−SS;
  current RDF + legacy both OS-only (SS gap). 2025 lumi already = 2.59933 nb⁻¹
  (CSV-verified); T_AA 2023 placeholder. RAA_plotting stale (Mac paths, legacy
  filenames, pT_bins_80 rebins). Committed all prior reco work (3 commits) first.

## Remaining Work
- **Year-combination normalization (NEEDS USER DECISION):** combined R_AA (and
  the existing combined crossx plots) naive-sum per-year `1/L_year`-weighted
  histos → absolute scale inflated by ~Σ(years)≈3. Fix consistently across
  crossx + R_AA (combine raw yields / total lumi). Until then R_AA absolute
  scale is preliminary (also 2023 T_AA placeholder).
- Replace reco placeholder with proper 3D pair ε_reco when Run 3 MC lands.
- Optional: fully retire `SingleBAnalysis/` (RDF now produces all it did + SS).
- Official 2024/2025 ⟨T_AA⟩ (currently 2023 placeholder).

## COMPLETION SUMMARY (2026-06-16) — doc CLOSED
Reco-corrected R_AA now runs from the RDF crossx, combined over PbPb 23+24+25 vs
pp24, with OS−SS combinatorial subtraction. Deliverables:
- SS signal-region histos added to RDF crossx (PbPb 3D, pp 2D); `/review-analysis-code` PASS.
- `RAA_plotting.cxx` modernized (cluster paths, RDF inputs case 6, OS−SS for
  pp+PbPb, combined years, 15-bin pT, index-wise ratio, mode-3 off-by-one fixed,
  segfault fixed); `/review-analysis-code` PASS (iter 2).
- R_AA plots vs pair pT / pair η / centrality in
  `dimuon_data/plots/single_b_analysis/RAA/`; `/review-plot` PASS (iter 4).
- 2025 lumi verified (2.59933 nb⁻¹, CSV-matched); 2023 T_AA placeholder kept.
Logs: `.claude/logs/review-analysis-code-20260616-013524-ss-crossx-histos-for-raa.md`,
`...-015446-raa-plotting-modernize.md`, `.claude/logs/review-plot-20260616-020322-raa-reco-corrected-plots.md`.

## Latest Stage
*(cleared — work complete; see COMPLETION SUMMARY. One open user decision:
year-combination normalization, see Remaining Work.)*
