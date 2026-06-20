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
Common notation (Run 2 HF-muon paper HION-2019-58 / arXiv:2109.00411, Eq. for
R_AA): **n_AA is the RAW yield (a count), not a per-event yield**; the per-event
normalization is carried by the explicit 1/N_evt factor.
```
R_AA(X) = 1/(⟨T_AA⟩ · N_evt) · (dn_AA/dX) / (dσ_pp/dX),   X ∈ {pair p_T, pair η}, per centrality
```
- **n_AA**: RAW OS−SS signal-region yield in PbPb (a count, efficiency-corrected
  per pair, then OS−SS subtracted bin-by-bin). NOT per-event-normalized.
- **N_evt**: number of sampled minimum-bias events in the centrality class =
  `f_ctr · σ_PbPb · L_PbPb`. In the code this is NOT counted separately — it is
  carried by the luminosity inside `crossx_factor` (see below), because
  `crossx_factor = 1/(f_ctr·σ_PbPb·⟨T_AA⟩·L_PbPb) = 1/(N_evt·⟨T_AA⟩)`. So one
  `weight_for_RAA` factor supplies BOTH 1/N_evt and 1/⟨T_AA⟩ at once.
- **dσ_pp/dX**: pp differential cross-section = (1/L_pp)-weighted OS−SS yield.
- **⟨T_AA⟩**: Glauber nuclear overlap per centrality. **2023 placeholder** for
  all years (roadmap Q2.3); folded into `weight_for_RAA`/crossx_factor already.
- Both pp and PbPb yields are efficiency-corrected (trigger × reco placeholder).
- The TAA-weighted crossx histogram IS the R_AA numerator
  `1/(⟨T_AA⟩·N_evt)·dn_AA/dX`; plot y-titles use the raw-yield symbol `n_AA`.

### 3. Step-by-step
a. **Yield & combinatorial background (§4a):** count OS pairs in the signal
   region (minv∈[1.08,2.9], pair p_T>8, −2.4≤q·η<2.2 per muon [one-sided, NO abs:
   forward single-muon acceptance is asymmetric in q·η], dR>0.05). Subtract
   the uncorrelated combinatorial component estimated from **SS** pairs passing
   the *same* cuts: `dN/dX = OS − SS`, bin by bin, BEFORE forming R_AA.
b. **Efficiency:** OS and SS both efficiency-corrected with the same per-muon
   weights (trigger × reco placeholder), so the subtraction is of corrected yields.
c. **Normalization:** PbPb via crossx_factor = 1/(σ_PbPb·f_ctr·⟨T_AA⟩·L_PbPb)
   (per centrality, in `weight_for_RAA`); this equals 1/(N_evt·⟨T_AA⟩), i.e. it
   supplies BOTH the 1/N_evt per-event normalization AND 1/⟨T_AA⟩. pp via 1/L_pp
   (`crossx_weight`).
d. **R_AA:** divide normalized PbPb by normalized pp per (X, centrality). Both
   1/N_evt and 1/⟨T_AA⟩ are already in the weight (via crossx_factor) — do NOT
   multiply by 1/⟨T_AA⟩ or 1/N_evt again.

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
- 2026-06-16 — Step 7 (year-combination fix) IMPLEMENTED. Confirmed correct vs
  HF R_AA note §4.1 Eq.3 (luminosity-weighted average). New header
  `Utilities/PbPbSampledLumi.h` (`PbPbMu4SampledLumiNb`: 23=1.17576
  [corrected 2026-06-19, GRL+b-hadron runs 461674+462964; was 1.02426], 24=0.85112
  [corrected 2026-06-19, GRL ≥489703; was 1.59663], 25=2.59933 nb⁻¹). Lumi-weighted combine `Σ(L_y·h_y)/ΣL` applied in: R_AA
  `HistRetrieve` (PbPb OS+SS 3D); crossx combined plotter
  `SingleBCrossxPlotterPbPbCombined::GetHistObject` (cross-section histos;
  `_counts` kept as simple sum); before/after stage plotter `getStage2D`.
  Reran R_AA + crossx + stage plots, all clean. **R_AA scale verified physical:**
  was ~1.5–2.5 mid-pT → now ~0.5–0.8 (÷~3 = ÷N_years in equal-σ limit; ΣL=5.22 nb⁻¹).
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
- ~~Year-combination normalization~~ **DONE 2026-06-16** (Step 7): luminosity-
  weighted average (HF R_AA note Eq.3) across crossx combined plotter + R_AA +
  stage plotter via `Utilities/PbPbSampledLumi.h`. `/review-*` PASS. R_AA scale
  now physical (~0.1–0.9).
- Official 2024/2025 ⟨T_AA⟩ (currently 2023 placeholder) → final R_AA scale.
- Citable σ_PbPb at 5.36 TeV (currently 7.8 b unvalidated guess) → R_AA scale.
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
*(CLOSED 2026-06-16 — task_06 complete incl. Step 7 year-combination fix; see
Step 7 DONE below. Record of the Step 7 plan retained for history.)*

**Step 7 (DONE): fix year-combination normalization** (user
approved; confirmed correct by the HF-muon R_AA note).

**Procedure (HF R_AA note HION-2019-58, §4.1 Eq.3):** "The central value is
extracted from the **luminosity-weighted average** of the fitted yields in 2015
and 2018 Pb+Pb data ... central value closer to 2018 due to its larger
luminosity." → combined = Σ(L_year·h_year)/ΣL_year. For our per-year histos
already weighted by crossx_factor (∝1/L_year), this equals the correct
ΣN/(f·σ·T_AA·ΣL). Weights = PbPb HLT_mu4 prescale-corrected sampled lumi:
2023=1.17576 (corrected 2026-06-19; was 1.02426), 2024=0.85112 (corrected 2026-06-19; was 1.59663), 2025=2.59933 nb⁻¹ (PbPbBaseClass crossx factors).

**Plan (Step 7):**
1. New single-source header `Utilities/PbPbSampledLumi.h` → `PbPbMu4SampledLumiNb(year)`.
2. R_AA `HistRetrieve`: lumi-weight the PbPb year combine (Σ L·h / ΣL); pp single year unchanged.
3. Crossx combined plotter `SingleBCrossxPlotterPbPbCombined::GetHistObject`:
   lumi-weighted average for cross-section histos; `_counts` histos = simple sum.
4. Before/after stage plotter `plot_crossx_reco_eff_stages.C` `getStage2D`: same
   lumi-weighting (consistency).
5. Rerun R_AA + crossx plots + before/after; `/review-analysis-code` + `/review-plot`.
6. Verify R_AA absolute scale drops by ~Σ(years) → physical (~0.5-0.8).

### Step 7 — DONE (2026-06-16): doc CLOSED
Year-combination normalization fixed to luminosity-weighted average (HF R_AA
note Eq.3) in all 3 combine sites (R_AA, crossx combined plotter, stage plotter)
via single-source `Utilities/PbPbSampledLumi.h`. `/review-analysis-code` PASS;
`/review-plot` PASS — R_AA absolute scale now physical (~0.1–0.9, was ~1.5–2.5).
Remaining caveats are physics inputs only: 2023 T_AA placeholder, σ_PbPb guess,
Run 2 reco placeholder — all flagged in `placeholder.md`.

*(Latest Stage cleared — task_06 complete.)*

---

## REOPENED 2026-06-19 — 2024 lumi correction (GRL ≥489703) + R_AA equation notation

**Trigger (user):** (1) PbPb 2024 lumi was computed with an old GRL that did not
exclude bad runs before 489703 → over-counted. Correct it and propagate to R_AA.
(2) The R_AA equation in this doc used non-common notation (treated dN as
per-event); fix to common notation (n_AA = raw yield, explicit 1/N_evt; cf. Run 2
HF-muon paper) and update the TAA-weighted crossx plot y-title accordingly.
(3)+(4) Explanations (delivered in chat): what enters R_AA; where TAA-weighted
SS/OS histos are filled.

### Finding (scope) — bad-run EVENTS are NOT in the numerator
`SkimCode/run_24hi/grid_sub.sh` shows the 2024 skim already used
`physics_HI2024_50ns.xml` with **runs ≥489703 only** ("pre-GRL runs removed";
datasets.tex event count 92.65M is post-GRL). So only the **luminosity CSV** was
over-counted (computed via lumicalc with a too-permissive GRL). ⇒ Pure
normalization fix, no re-skim.

### Corrected value
Drop the 19 pre-489703 runs from `lumitable_pbpb_24_HLT_mu4.csv`. New Total
(Prescale Corrected) = **851.118101 µb⁻¹ = 0.85112 nb⁻¹** (was 1596.63 µb⁻¹ =
1.59663). Other Total cols: Good 9455→4401, LDelivered→877.008663,
LRecorded→851.118721, LAr Corrected→851.118101, Live 96.60→97.05, LAr→100.00.

### Net physics effect (combined R_AA / TAA-weighted crossx)
L cancels in `L_y·h_y`, so the only change is the year-combine denominator ΣL:
5.22022 → 4.47471 nb⁻¹ ⇒ **combined R_AA and TAA-weighted crossx × 1.1666
(+16.7%)**. (Per-2024-year crossx, if ever viewed alone, scales ×1.876 = 1/0.5328.)

### Plan (this rework)
1. **CSV** `IntNotes/data/luminosity/pbpb_2024/lumitable_pbpb_24_HLT_mu4.csv`:
   delete runs <489703, recompute Total row. (mu4_mu4noL1 is reference-only /
   not maintained — flag to user, update only if requested.)
2. **Code (coupled, /review-analysis-code):** `PbPbBaseClass.h`
   make_crossx_factors_pbpb_2024 (6× 1.59663→0.85112 + comment) AND
   `Utilities/PbPbSampledLumi.h` (24: 1.59663→0.85112 + comment). Both MUST
   change together (combine cancellation). Refill required.
3. **Docs:** `IntNotes/analysis_metadata.md`, `IntNotes/tex/datasets.tex`,
   `IntNotes/data/luminosity/README.md`, `Analysis/docs/placeholder.md`,
   `analysis_roadmap_2026_06.md` — 2024 lumi 1.59663→0.85112; note GRL ≥489703.
4. **Notation (/review-plot covers relabel):** equation in this doc fixed (DONE
   above, §2/§3c/§4); y-title in `plot_single_b_crossx_pbpb.cxx` d²N→d²n_{AA},
   dN→dn_{AA} (keep 1/(⟨T_AA⟩ N_evt) prefactor).
5. **Rerun:** recompile PbPb RDF; rerun 2024 PbPb crossx (refill); rerun combined
   crossx plotter + RAA_plotting. /review-plot on regenerated plots. Verify
   combined R_AA ×~1.17.

### Progress (2026-06-19) — COMPLETE
- Step 1 CSV: 19 pre-489703 runs deleted from `lumitable_pbpb_24_HLT_mu4.csv`;
  Total Prescale Corrected = 851.118101 µb⁻¹ = 0.85112 nb⁻¹. (mu4_mu4noL1
  reference-only — NOT touched; flagged to user.)
- Step 2 code: PbPbBaseClass.h make_crossx_factors_pbpb_2024 (6× 0.85112) +
  PbPbSampledLumi.h case 24. `/review-analysis-code` PASS iter 1 (log
  review-analysis-code-20260619-010258-pbpb-2024-lumi-grl-correction.md).
  ACLiC clean; 2024 crossx refilled (run_crossx_hist_filling_pbpb24.sh).
- Step 3 docs: analysis_metadata.md (+2024 GRL subsection), datasets.tex,
  luminosity/README.md, placeholder.md, roadmap, this doc — all 0.85112.
- Step 4 notation: equation §2/§3c/§4 rewritten (n_AA raw yield, explicit
  1/N_evt). Plot y-titles d²N→d²n_{AA}, dN→dn_{AA} (counts keep N_events).
- Replot: combined PbPb crossx (pt_120 + pt_150) + R_AA regenerated.
  `/review-plot` PASS iter 1 (log
  review-plot-20260619-011739-pbpb-2024-lumi-raa-relabel.md). R_AA physical
  (~0.1–1.2), centrality ordering correct, +16.7% scale verified.

### Signal-cut consistency fix (2026-06-20) — q·η, one-sided
User flagged the `|q·η|<2.2` abs-value mislabel. The cut is per-muon **q·η<2.2,
ONE-SIDED** (−2.4≤q·η<2.2): forward single-muon acceptance is asymmetric in q·η
(toroid bends in the polar plane; charge×polarity sets the side). Two fixes
(`/review-analysis-code` PASS iter 1, log
review-analysis-code-20260619-235854-qeta-permuon-cut-fix.md):
- **FIX A (label only):** `CommonEffcyConfig.h:15,41,66` and §3a above
  `|q×eta|<2.2` → `−2.4≤q·η<2.2`. Binning arrays unchanged (already −2.4…2.2).
- **FIX B (physics):** fullsim/overlay signal cut used PAIR η (`pair_eta<2.2`),
  inconsistent with the per-muon data signal region. Changed to per-muon q·η on
  BOTH muons (truth: `m1.truth_charge*m1.truth_eta<2.2 && m2...`; reco:
  `m1.charge*m1.eta<2.2 && m2...`) in `RDFBasedHistFillingPythiaFullsim.cxx`
  (129/131), `…PythiaFullsimOverlay.cxx` (69/71), `…PowhegFullsim.cxx` (185/186).
  ACLiC-clean; samples NOT rerun (placeholder reco-eff; matters for future real
  3D pair ε_reco). Legacy `SingleBAnalysisBase.cxx:25` /
  `FillHistogramsCrossx_PP_clean.cxx:8` flagged (pair_eta), left (retired).

### Q3/Q4 (explanations) — delivered in chat 2026-06-19
Q3 (what enters R_AA) + Q4 (where TAA-weighted SS/OS histos are filled, with
code line refs) answered. Key line refs (RDFBasedHistFillingPbPb.cxx): weight
940/954/955/973; OS 3D 1006-1009; SS 1027-1043; per-ctr 2D 1061-1145.

### Latest Stage
*(cleared — 2024 lumi correction + R_AA notation rework COMPLETE; both review
loops PASS. Changes uncommitted, awaiting user go-ahead to commit.)*

---

## REOPENED 2026-06-19 (2) — 2023 lumi correction (old GRL) + subtract two b-hadron runs

**Trigger (user):** PbPb 2023 lumi table was computed with an old (wrong) GRL.
User has replaced `IntNotes/data/luminosity/pbpb_2023/lumitable_pbpb_23_HLT_mu4.csv`
with the correct-GRL table and deleted the obsolete 2023 `mu4_mu4noL1` CSV.
Tasks: (1) update the 2023 lumi README + other relevant docs; (2) for the
luminosity used in R_AA, subtract the **two b-hadron runs** — previously only ONE
(462964) was in the old table's run range and excluded; the new (larger) GRL table
now contains BOTH, so both must be subtracted.

### Corrected value (§2/§3c normalization, L_PbPb for 2023)
- New GRL total (Prescale Corrected, summed over runs) = **1183.650457 µb⁻¹ =
  1.18365 nb⁻¹** (old GRL gave 1029.52 µb⁻¹).
- Two b-hadron runs to subtract: **461674** (2.623332 µb⁻¹) + **462964**
  (5.262407 µb⁻¹) = 7.885739 µb⁻¹.
- **R_AA luminosity (2023) = 1183.650457 − 2.623332 − 5.262407 = 1175.764718
  µb⁻¹ = 1.17576 nb⁻¹** (was 1.02426; old value = old-GRL total 1029.52 − 5.262
  for run 462964, with 461674 below the old table range).

### Finding (scope) — b-hadron-run EVENTS are EXCLUDED from the numerator (consistent)
The two b-hadron runs are dropped from the analysis at NTuple processing:
`PbPbEventSelConfig.h::PbPbBadRuns(23) = {461674, 462964}`, consumed by
`PbPbExtras::PassEventSel()` (PbPbExtras.c:206–209: `if RunNumber ∈ bad_runs
return false`), so no event from either run enters the muon_pairs trees that feed
the crossx RDF. (The 2023 skim GRL v120 `HeavyIon_All_Good`,
`SkimCode/run_23hi/grid_sub.sh`, does contain these runs, but the offline
event-selection bad-run cut removes them.) ⇒ Excluding their luminosity from the
denominator is the **physically consistent** matching move — both numerator
(events) and denominator (lumi) drop these runs. The old 1.02426 already excluded
462964's lumi; the only gap was 461674 (absent from the old table's run range but
present in the new GRL table), now also subtracted. No numerator change needed
(events already excluded); this is the denominator catching up to the event cut.

### Net physics effect (combined R_AA / TAA-weighted crossx)
L cancels in `L_y·h_y` (h_y ∝ 1/L_y), so only the year-combine denominator ΣL
changes: 4.47471 → 4.62621 nb⁻¹ ⇒ **combined R_AA and TAA-weighted crossx ×
0.96725 (−3.27%)**. (Per-2023-year crossx alone, if ever viewed: × 0.87114 =
1.02426/1.17576, −12.89%.)

### Plan (this rework)
1. **CSV** (user already done): correct-GRL `lumitable_pbpb_23_HLT_mu4.csv` in
   place; obsolete `mu4_mu4noL1` deleted. Verify total + b-hadron run values (DONE
   above). No re-skim (events already include these runs; lumi-only per user).
2. **Code (coupled, /review-analysis-code):** `PbPbBaseClass.h`
   make_crossx_factors_pbpb_2023 (6× 1.02426→1.17576 + comment) AND
   `Utilities/PbPbSampledLumi.h` (case 23: 1.02426→1.17576 + comment). Both MUST
   change together (combine cancellation). Refill required. Quote §2/§3c.
3. **Docs:** `IntNotes/data/luminosity/README.md` (2023 row + sanity anchor),
   `IntNotes/analysis_metadata.md` (2023 lumi + GRL/b-hadron note, drop
   mu4_mu4noL1 2023), `IntNotes/tex/datasets.tex`, `Analysis/docs/placeholder.md`,
   `analysis_roadmap_2026_06.md` (Q2.1), `analysis_status_summary.md` — 2023 lumi
   1.02426→1.17576.
4. **Rerun:** recompile PbPb RDF; rerun 2023 PbPb crossx (refill); rerun combined
   crossx plotter + RAA_plotting. /review-plot on regenerated plots. Verify
   combined R_AA × ~0.967.

### Progress (2026-06-19) — Steps 1–2 DONE, 3–4 in progress
- Step 1 CSV: correct-GRL `lumitable_pbpb_23_HLT_mu4.csv` in place (user); obsolete
  2023 `mu4_mu4noL1` CSV deleted (user). Verified: Prescale-Corrected data-row sum
  = 1183.650457 µb⁻¹; b-hadron runs 461674 (2.623332) + 462964 (5.262407) µb⁻¹.
- Step 2 code: PbPbBaseClass.h make_crossx_factors_pbpb_2023 (6× 1.02426→1.17576 +
  comment) + PbPbSampledLumi.h case 23 (1.17576 + comment). `/review-analysis-code`
  PASS iter 1 (log review-analysis-code-20260619-230240-pbpb-2023-lumi-grl-bhadron.md):
  1.17576 re-derived from CSV MATCH; coupling invariant held; 2024/2025 untouched;
  ACLiC clean (PbPbMu4SampledLumiNb 23/24/25 = 1.17576/0.85112/2.59933).
- Step 3 refill: `run_crossx_hist_filling_pbpb23.sh` rc=0 →
  `pbpb_2023/histograms_real_pairs_pbpb_2023_single_mu4_no_trg_plots_nominal.root`
  regenerated (per-2023 crossx × 0.87114 vs pre-correction; combined cancels).
- Numerator consistency CONFIRMED (Finding above corrected): both runs excluded at
  NTuple processing via PbPbExtras::PassEventSel/PbPbBadRuns(23) — events already
  out of the numerator, so subtracting their lumi is consistent.

- Step 4 docs: analysis_metadata.md (2023 GRL-correction + b-hadron-run-subtraction
  section rewritten; lumi row 1183.65 µb⁻¹/prescale 1.4312; code-sync status;
  prescale note), luminosity/README.md (dropped 2023 mu4_mu4noL1 row; physical
  total 1.18 nb⁻¹; sanity anchor Σ=4.63 nb⁻¹), datasets.tex (lumi table 1.17576,
  run-exclusion prose, caption ≃1.43), roadmap Q2.1 + status ledger (×2 entries),
  status_summary Latest Update. placeholder.md: no 2023-lumi mention → no change.
- Step 5 replot: combined PbPb crossx (TAA_weighted 30 + counts 36 PNGs) +
  R_AA (3 PNGs, vs pT/η/ctr) regenerated, rc=0. `/review-plot` PASS iter 1 (log
  review-plot-20260619-231947-pbpb-2023-lumi-crossx-raa.md): ΣL 4.47471→4.62621,
  combined scale 0.96725 verified; R_AA physical (~0.1–1.2), centrality-ordered;
  crossx combined-only; OS−SS applied; refilled 2023 op−ss=14490.3 (SS/OS≈28%).

### COMPLETION SUMMARY (2026-06-19 (2)) — 2023 lumi rework DONE
PbPb 2023 sampled luminosity corrected for a wrong/old GRL and the two b-hadron
runs (461674 + 462964) subtracted from the R_AA luminosity: **1.02426 → 1.17576
nb⁻¹** (correct-GRL total 1.18365 nb⁻¹ − 7.886 µb⁻¹ for the two runs, both already
excluded at event level via PbPbBadRuns/PassEventSel → numerator+denominator
consistent). Coupled change in PbPbBaseClass.h (6 factors) + PbPbSampledLumi.h
(case 23), `/review-analysis-code` PASS; 2023 crossx refilled; combined crossx +
R_AA replotted, `/review-plot` PASS. Net: combined R_AA / TAA-weighted crossx
**× 0.96725 (−3.3%)** (ΣL 4.47471→4.62621 nb⁻¹). Obsolete 2023 mu4_mu4noL1 CSV
deleted (user). Docs synced (metadata, lumi README, datasets.tex, roadmap, status).
Logs: review-analysis-code-20260619-230240-pbpb-2023-lumi-grl-bhadron.md,
review-plot-20260619-231947-pbpb-2023-lumi-crossx-raa.md.

### Latest Stage
*(cleared — 2023 lumi GRL correction + two-b-hadron-run subtraction COMPLETE; both
review loops PASS. Changes uncommitted, awaiting user go-ahead to commit — same
batch as the 2024 lumi correction + R_AA notation rework above.)*
