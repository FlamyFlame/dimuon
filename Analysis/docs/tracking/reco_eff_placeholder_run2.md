# Reco-efficiency placeholder from Run 2 internal notes

**Mode:** Implementation. **Created:** 2026-06-15.
**Reviewer rules:** code steps → `/review-analysis-code`; plot steps → `/review-plot`.
**Sibling docs:** `analysis_roadmap_2026_06.md` (§Q4, chain step [4]/[5], task_05),
`placeholder.md` (registry item 3), `analysis_status_summary.md` (crossx status).

---

## Objective

The Run 3 reconstruction-efficiency MC (Pythia fullsim pp24 + Pythia fullsim
HIJING overlay) is unavailable for ≥2 months. Install a **temporary Run 2
single-muon reco-efficiency placeholder** so the crossx (and later R_AA)
chain runs end-to-end with a reco-eff correction applied, while keeping the
correction clearly labelled as a placeholder to be replaced by the proper
Run 3 **pair** efficiency ε_reco(pair pT, pair η, dR).

Deliverables:
0. Choose pp + PbPb placeholders (DONE — see Physics Procedure §3).
1. Digitize the Run 2 single-muon reco efficiencies into a ROOT file;
   reproduce the F.2 plots; `/review-plot` cross-check vs the originals.
2. RDF correction-stage **enum** (unfolding / reco-eff / trig-eff) in pp +
   PbPb, applied sequentially; save histograms: raw → unfolded →
   unfolded+reco → unfolded+reco+trig.
3. Back up current pp & PbPb crossx plots (incl. before/after eff plots).
4. Apply ε₁·ε₂ reco placeholder weight in the crossx RDF; rerun crossx
   (pp + PbPb, skip ntuple). Before/after eff comparison = 3 lines:
   raw, reco-only (placeholder), reco+trig.
5. Document the choice in `placeholder.md`; update roadmap + status docs.

---

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation
The differential cross-section divides the raw signal yield by the net
per-pair efficiency to recover the true production rate. Without a
reconstruction-efficiency correction the cross-section is biased low
(muons that were produced but not reconstructed are missing), most
severely at low muon pT where reco efficiency turns on. Until Run 3 MC
exists we substitute Run 2 efficiencies as a stopgap.

### 2. Top-level equation (per opposite-sign signal pair)
The net per-pair correction weight applied to the yield is
```
w⁻¹ = ε_trig^pair(p_a, p_b) · ε_reco(p_a) · ε_reco(p_b)
```
- `ε_trig^pair` — per-pair trigger efficiency, **already applied** in the
  crossx RDF: PbPb single-mu4 inclusion–exclusion ε₁^nc+ε₂^nc−ε₁^nc·ε₂^nc;
  pp24 2mu4 ε₁^nc·ε₂^nc (`RDFBasedHistFillingPbPb.cxx:937`, `PP.cxx:331`).
- `ε_reco(p)` — **single-muon** reconstruction efficiency as a function of
  (pT, q·η[, centrality]). The pair reco efficiency is taken as the
  **product** ε_reco(p_a)·ε_reco(p_b) — this is the Run 2 dimuon-note
  treatment (Eq. eq:net_eff_correction, ATL-COM-PHYS-2021-1094 §Corrections).
- p_a, p_b are the two muons; q·η = charge × η per muon.

**This product is a PLACEHOLDER proxy and is physically wrong for our
signal** (the two signal muons are nearby; pair reconstruction is
correlated, so it does NOT factorize into single-muon products — see
`docs/tracking/hijing_overlay_reco_effcy_investigation.md` and
`atlas_inner_detector_tracking` KB on close-track merged clusters). The
proper Run 3 correction is the 3D pair efficiency ε_reco(pair pT, pair η,
dR) of task_05 / roadmap §Q4, applied once Pythia fullsim + overlay land.

### 3. Single-muon ε_reco sources (Run 2 placeholders)

**PbPb — ATL-COM-PHYS-2021-1094 (Run 2 dimuon note), Appendix F.2.**
Single-muon ε_reco for **Medium** muons vs pT^truth, in 9 q·η slices, per
centrality interval, from STARlight+HIJING-overlay MC at √s_NN = 5.02 TeV.
- q·η slices (match `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap_run2`):
  [-2.4,-2.0],[-2.0,-1.5],[-1.5,-1.0],[-1.0,-0.5],[-0.5,0.5],[0.5,1.0],
  [1.0,1.5],[1.5,2.0],[2.0,2.4].
- centrality intervals (7): 0-10,10-20,20-30,30-40,40-50,50-60,60-80%.
- ε rises from ~0.5–0.65 at pT≈4 GeV to a ~0.8–0.9 plateau by pT≳8 GeV;
  small systematic increase central→peripheral.
- WP = **Medium** (F.2). F.1 (Tight) is NOT used — the user pointed to F.2.

**pp — HION-2019-58 / arXiv:2109.00411 (Run 2 HF-muon R_AA note), Fig. 31.**
Data-driven single-muon ε_reco^pp for **Medium** muons vs pT (MuonCP, 2017
13 TeV, applied to 5.02 TeV), in two |η| ranges:
- barrel 0.10<|η|<1.05: ε ≈ 0.77 (pT 4–5) → ~0.93 (5–6) → ~0.96 (6–9) →
  ~0.98 plateau (pT≳10).
- endcap 1.30<|η|<2.10: ε ≈ 0.60 (4–5) → ~0.85 (5–6) → ~0.93 (6–8) →
  ~0.96–0.97 plateau.
- pp has **no centrality**. This is the SAME source the Run 2 dimuon note
  cites for its own pp efficiencies (§Corrections, "pp efficiencies taken
  from HION-2019-58"), hence the preferred pp placeholder.

**pp fallback (NOT used):** if no pp source had been found, the most
peripheral (60–80%) PbPb F.2 curve would be the best alternative (lowest
hadronic activity → closest to pp). A pp source WAS found, so this fallback
is unused; recorded for provenance.

### 4. Centrality mapping (F.2 7 intervals → analysis 6 bins)
Analysis crossx centrality bins (`FindCtrSuffix`): 0-5,5-10,10-20,20-30,
30-50,50-80. F.2 intervals are mapped by containment / averaging:
- 0-5, 5-10 → F.2 0-10
- 10-20 → 10-20; 20-30 → 20-30
- 30-50 → mean(F.2 30-40, 40-50)
- 50-80 → mean(F.2 50-60, 60-80)
The lookup picks the F.2 interval whose range contains the muon's event
centrality (using the F.2 edges directly), so no pre-averaging is needed in
the stored file; the file stores all 7 F.2 intervals.

### 5. Negative constraints (what the code must NOT do)
- MUST NOT treat ε₁·ε₂ as the final/correct pair efficiency — it is a
  flagged placeholder. Every code site carries a `// PLACEHOLDER` remark
  pointing to task_05 / the proper 3D pair efficiency.
- MUST NOT reuse the trigger-efficiency fit files for reco — reco eff is a
  different quantity (offline reco vs HLT match) with a different source.
- MUST NOT apply any MuonCP scale factor on top (roadmap Q2.7: omitted).
- MUST NOT silently overwrite the existing trigger-corrected histograms;
  the stage histograms are additive (new outputs), and with all reco ε set
  to 1 the trig+reco histogram must equal the existing trig-only result.
- "unfolded" stage is currently an **identity** (no unfolding input yet,
  roadmap Q4) — it is a framework hook, equal to raw until unfolding lands;
  must be labelled as such.

---

## Context
- Run 2 method text: `IntNotesRun2DimuonReference/tex/Corrections.tex`
  (§reco_eff) + `tex/Appendix_RecoEffs.tex` (F.1 Tight, F.2 Medium).
- Run 2 dimuon note PDF: `IntNotesRun2DimuonReference/ATL-COM-PHYS-2021-1094.pdf`
  (F.2 = pp. 97–103). HF R_AA note PDF: `.claude/kb/analysis/Run2 HF muons
  nuclear modification factor internal notes.pdf` (Fig. 31 = pp. 41–42).
- Existing per-muon lookup pattern: `RDFBasedHistFillingData.cxx`
  `EvaluateSingleMuonEffcy*`, `OpenEffcyPtFitFile`, `FindCtrSuffix`.
- Crossx weight sites: PbPb `RDFBasedHistFillingPbPb.cxx:913+`
  (`weight_for_RAA_trig_corr`), pp `RDFBasedHistFillingPP.cxx:348+`
  (`crossx_weight_trig_corr`).
- q·η + pair-eta binning: `RDFBasedHistFilling/CommonEffcyConfig.h`.

## Scope
In: digitized Run 2 reco-eff ROOT file + reproduction plots; RDF stage-enum
+ reco placeholder weight (pp + PbPb crossx path only); crossx rerun (skip
ntuple); before/after 3-line plots; docs. Out: R_AA (task_06), unfolding,
the proper 3D pair efficiency (task_05 proper), systematics.

## Design Decisions
- **Stage histograms added only for the primary `pair_pt × pair_eta`
  differential** (not every observable) — keeps output size bounded while
  giving the pT and η differential cross-sections at each correction stage,
  which is what the before/after comparison needs. Other observables keep the
  existing single (trig-only) weight. Rationale: minimal scope; the staged set
  is for impact visualization, not a full re-derivation of every histogram.
- **Backward compatibility preserved:** existing `..._w_signal_cuts`
  (trig-only), `..._no_trig_corr`, `..._counts` histograms unchanged; the
  4 stage histograms are ADDED with `_corr_{raw,unfolded,unfolded_reco,
  unfolded_reco_trig}` suffixes. Invariant (Physics Procedure §5): with
  w_reco=1 and w_unfold=1, `_corr_unfolded_reco_trig` == existing trig-only.
- **w_reco lookup reuses no trig machinery** (separate `s_reco_eff_ph_map`
  TGraph store + `OpenRecoEffPlaceholderFile`/`EvaluateSingleMuonRecoEffPlaceholder`)
  per Physics Procedure §5 negative constraint.
- **1/ε guard:** pair efficiency floored at 0.05 (w_reco ≤ 20); lookup failure
  → w_reco=1 (no correction, never drops the pair).
- **pp centrality sentinel:** pp calls pass centrality=-1 → barrel/endcap by
  |q·η|<1.05; PbPb maps avg_centrality onto the F.2 7 intervals.
- **enum/stage table shared** via `CorrectionStages.h`; same cumulative column
  names (`cw_*`) in PP and PbPb so the table drives both.

## Implementation Plan
1. Tracking doc + Physics Procedure. — per §all. (this file)
2. Digitize F.2 (PbPb) + Fig.31 (pp) → ROOT file of ε_reco vs pT per
   (q·η, centrality / η-region). Per §3. Builder is a small ROOT macro.
3. Reproduce single-muon eff plots (F.2 layout) from the ROOT file →
   `/review-plot` vs the original PDF panels. Per §3.
4. RDF: correction-stage enum + sequential weights + stage histograms,
   pp + PbPb. Per §2/§5. → `/review-analysis-code` (quote §2,§5).
5. RDF: reco placeholder lookup `Evaluate...RecoEffcyPlaceholder` +
   ε₁·ε₂ weight + 1/ε clamp. Per §2,§4,§5. → `/review-analysis-code`.
6. Back up current pp & PbPb crossx plots.
7. Rerun crossx pp + PbPb (SKIP ntuple); sanity-check outputs.
8. Before/after eff plots (3 lines: raw, reco-only, reco+trig) →
   `/review-plot`.
9. Update `placeholder.md`, `docs/`, roadmap ledger, memory.

## Progress Log
- 2026-06-15 — Step 1: doc created. Confirmed (a) Run 2 dimuon note uses
  ε_reco(p_a)·ε_reco(p_b) single-muon product (Corrections.tex), justifying
  the ε₁·ε₂ proxy; (b) F.2 = Medium muons, 9 q·η slices = coarse_incl_gap_run2,
  7 centrality intervals; (c) pp placeholder = HF R_AA Fig.31 Medium
  ε_reco^pp (barrel/endcap), the source the dimuon note itself cites for pp;
  (d) analysis crossx centrality bins = {0-5,5-10,10-20,20-30,30-50,50-80}.

## Results & Observations
- **Digitized values** (anchor pT = 4,4.5,5,6,7,8,10,13,16,19 GeV). Base
  (0–10%) per q·η slice + additive centrality bump {0,.01,.02,.03,.035,.04,.05}
  for {0-10..60-80}. Hardcoded in
  `plotting_codes/reco_effcy/build_run2_reco_eff_placeholder.C`.
- **Output ROOT:** `EfficiencyCorrs/EffFiles/run2_reco_eff_placeholder.root` —
  63 PbPb TGraphs `gr_reco_eff_medium_pbpb_ctr{lo}_{hi}_q_eta_{suffix}`
  (7 ctr × 9 q·η) + `gr_reco_eff_medium_pp_barrel`, `..._pp_endcap`.
- **Reproduction plots:** `plots/reco_effcy/placeholder/` (7 PbPb 3×3 + 1 pp).
  Visual match to F.2: flat forward/backward slices ~0.78–0.85; turn-on for
  |q·η|<1.5 from ~0.45–0.62 at pT=4 to ~0.85–0.95 plateau.

## Remaining Work
- **Pending user decision (not started):** promote the reco-inclusive
  `_corr_unfolded_reco_trig` stage to the **nominal** crossx + R_AA plots
  (`plot_single_b_crossx_{pp,pbpb}.cxx`, `RAA_plotting.cxx`). NOT done on
  purpose — this placeholder is known-poor; surfaced to the user. Standard
  crossx/R_AA plots currently still read the trigger-only histograms.
- Replace the whole placeholder with the proper 3D pair ε_reco(pair pT, pair η,
  dR) once Pythia fullsim (pp24) + HIJING overlay (r17662) land (task_05 proper).
- `plot_crossx_trig_corr_sanity.C` has a stale `DatasetTriggerMap::Get` call
  (correct API is `GetTrigger`); the new `plot_crossx_reco_eff_stages.C` avoids
  it by hardcoding trig suffixes. Fix the old macro if it is ever re-run.

## REOPENED 2026-06-16 — Follow-ups Q1 (PbPb differential crossx) + Q2 (MC-data auto-update)

### Q1 — PbPb differential cross-section plots (not T_AA-weighted)
**Finding:** the PbPb crossx `TAA_weighted/` plots are filled with
`weight_for_RAA_trig_corr` (= weight·crossx_factor·w_reco·w_trig, crossx_factor
= 1/(f·σ_PbPb·T_AA·L)) — i.e. the T_AA-weighted R_AA-input yield dN_AA/T_AA — but
MIS-LABELLED "d²σ/dp_Tdη [pb GeV⁻¹]". There is no genuine differential
cross-section (dσ/dp_T = (1/L)·dN/dp_T) in the PbPb crossx.
**Plan:** add a PbPb differential-cross-section weight in the RDF:
`weight_for_xsec_trig_corr = weight·(1/L_year)·w_reco·w_trig` (1/L via
`PbPbMu4SampledLumiNb`, nb units), filled per-centrality (+ no-corr variant
mirroring pp's `crossx_weight`). Plot a `differential_crossx/` set (dσ/dp_T
[nb GeV⁻¹], lumi-combined). Keep the T_AA-weighted 3D (R_AA input) but relabel
its plots honestly (dN_{AA}/T_{AA}) or drop the misleading dir. pp already plots
true dσ (crossx_weight = weight/L_pp) — no pp change. → /review-analysis-code + /review-plot.

### Q2 — MC-data comparison must reflect reco (and any future) efficiency/det-resp corrections
**Finding:** `plot_mc_data_compr.cxx` reads the **generic gapcut** data histos
`h_{kin}_{op,ss}_wgapcut` from the pp crossx output. Those are filled by
`FillHistogramsGeneric` with `generic_weight_col="w_trig"` — TRIGGER ONLY, NOT
reco-corrected. So the MC-data comparison does NOT currently reflect the reco
placeholder, and it is STALE (plots 06-10 vs pp crossx 06-16).
**Plan:** (a) fold w_reco into the pp generic analysis weight (in
`RDFBasedHistFillingPP::FillHistogramsGeneric`, the `!trigger_effcy_calc` branch:
define w_reco, set `generic_weight_col` to w_reco·w_trig) so the gapcut data
histos are reco+trig corrected; (b) rerun pp crossx RDF (regenerates generic
histos) + `plot_mc_data_compr.cxx`; (c) document the INVARIANT: MC-data
comparison MUST be regenerated after any pp efficiency/det-response/unfolding
change (the pp crossx pipeline stage 8 does this; direct RDF reruns must also run
it). → /review-analysis-code + /review-plot.
**Auto-update mechanism:** pp crossx pipeline `pipeline_pp_crossx.sh` stage 8 already
runs the comparison; the invariant note + (optional) wiring ensures it isn't skipped.

### Q1 + Q2 — DONE & REVIEWED (2026-06-16)
- **Q1:** PbPb genuine differential cross-section `weight_for_dsigma_trig_corr`
  = weight·(1/L_year)·w_reco·w_trig added to `RDFBasedHistFillingPbPb.cxx`
  (`PbPbSampledLumi.h`); per-centrality `*_dsigma` 2D (pair_eta/minv/dr) + 3D
  (dr×eta×pt) histos. `plot_single_b_crossx_pbpb.cxx`: `TAA_weighted/` plot block
  → `differential_crossx/` reading the `_dsigma` histos, labeled "dσ/dp_T
  [nb GeV⁻¹]". T_AA-weighted 3D still produced (R_AA input). Verified dσ/T_AA-weighted
  = f·σ·T_AA = 10.2 (ctr0_5, exact). `/review-analysis-code` + `/review-plot` PASS.
- **Q2:** reco folded into the pp generic analysis weight —
  `RDFBasedHistFillingPP::FillHistogramsGeneric` (`!trigger_effcy_calc` branch)
  now defines w_reco and sets `generic_weight_col="w_reco_trig"` (was "w_trig"),
  so the gapcut histos `h_{kin}_{op,ss}_wgapcut` read by the MC-data comparison are
  reco+trig corrected. Collision fix: reco Defines moved into the
  `generic_weight_col.empty()` guard in the crossx OS+SS blocks. Reran pp crossx
  RDF + `plot_mc_data_compr.cxx` (updated). `/review-*` PASS.
- **INVARIANT (Q2c):** the MC-data comparison (`plot_mc_data_compr.cxx`) MUST be
  regenerated after ANY change to the pp crossx efficiency / detector-response /
  unfolding corrections — because its data histos carry those corrections via
  `generic_weight_col`. This holds for: (1) this reco placeholder, (2) the future
  real efficiency corrections, (3) unfolding. The pp crossx pipeline
  `pipeline_pp_crossx.sh` stage 8 runs it automatically; when rerunning the pp
  RDF directly (`run_crossx_hist_filling_pp24.sh`), also rerun the comparison.
- **pt_150 variant — FIXED (2026-06-16):** added `_dsigma` pt_150 histos
  (`h2d_op_crossx_dsigma_vs_pair_eta_vs_pt_150_<ctr>`,
  `h3d_crossx_dr_vs_pair_eta_vs_pt_150_dsigma_<ctr>`) and switched the plotter's
  pt_150 block from `TAA_weighted/` to `differential_crossx/` (nb/GeV). Reran
  PbPb crossx + plotter (use_pt_bins_150=true); renders correctly. Removed the
  two now-orphaned live `TAA_weighted/` dirs (pt_120 + pt_150); backups kept.
  Mechanical mirror of the reviewed pt_120 dσ pattern (compile+render verified).
- **Remaining (follow-up):** optional L==0 guard on `dsigma_lumi_factor`
  (unreachable); proper 3D pair ε_reco when Run 3 MC lands.

## COMPLETION SUMMARY (2026-06-16)
[superseded by the Q1/Q2 follow-ups above; original placeholder work complete.]
All requested work done; tracking doc CLOSED (removed from Active list, file
kept). Deliverables:
- **Placeholder choices:** PbPb = dimuon note F.2 (Medium single-μ ε_reco per
  centrality); pp = HF R_AA note Fig.31 (Medium ε_reco^pp barrel/endcap). pp
  fallback (peripheral PbPb) recorded as unused. Memory:
  `project_pp_reco_eff_placeholder`.
- **ROOT file + reproduction plots** (`/review-plot` PASS).
- **Correction-stage enum framework** (`CorrectionStages.h`) + ε₁·ε₂ reco weight
  in PP+PbPb crossx RDF (`/review-analysis-code` PASS); 4 stage histograms.
- **Backups** of pre-rerun crossx plots + histogram ROOTs (`*_backup_20260615`).
- **Crossx reran** pp + 3 PbPb years (skip ntuple); **before/after 3-line plots**
  (`/review-plot` PASS).
- **Docs:** `placeholder.md` item 3, roadmap ledger, this doc; PLACEHOLDER
  remarks at every code site.

## Latest Stage
*(cleared — Step 10 propagation complete; doc RE-CLOSED 2026-06-16)*

### Step 10 — DONE (2026-06-16): reco placeholder promoted to NOMINAL
- Backup: `crossx_hist_backup_20260616_pre_reco_nominal/` (pre-reco trig-only
  histos + `RAA_plotting.cxx.bak`); plot dirs `*_backup_20260616_pre_reco_nominal`.
- `w_reco` folded into nominal weights: PbPb `weight_for_RAA_trig_corr` =
  `weight_for_RAA·w_reco·w_trig`; pp `crossx_weight_trig_corr` =
  `crossx_weight·w_reco·w_trig`. ACLiC clean. `/review-analysis-code` PASS
  (`.claude/logs/review-analysis-code-20260616-010859-reco-fold-nominal-weight.md`):
  algebra == cw_unfolded_reco_trig, no double-count, ordering valid.
- Crossx RDF reran (pp + 3 PbPb, rc=0). **Verified** nominal == reco_trig stage:
  pp nominal 5739 (was 4367 trig-only, ×1.314); PbPb23 R_AA 3D input ×1.665;
  ctr10_20 nominal == stage.
- Nominal crossx plots reran: `plots/single_b_analysis/{pp24,pbpb_23_24_25_combined}`.
  Self-verified render OK; content == already-`/review-plot`-validated reco_trig
  stage (plotter code unchanged → no redundant re-review).
- **R_AA:** `RAA_plotting.cxx` stale legacy (Mac `base_dir`, `*_single_b_ana_hists*`
  filenames, reads non-produced `h3d_ss_crossx_...`) → NOT runnable; modernization
  = task_06. Reco placeholder IS in the R_AA input histogram now.

### (history) original reopen note
**REOPENED 2026-06-16 (user request):** propagate the reco placeholder to the
**nominal** crossx + R_AA (Run 3 MC won't land for 2–3 months; downstream steps
must proceed on reco-corrected nominal).

**Plan (Step 10):**
1. Fresh backup (today) of current nominal crossx histos + plots + `RAA_plotting.cxx`
   before overwriting.
2. **Propagation mechanism:** fold `w_reco` into the nominal corrected-weight
   columns — PbPb `weight_for_RAA_trig_corr` = `weight_for_RAA * w_reco * w_trig`;
   pp `crossx_weight_trig_corr` = `crossx_weight * w_reco * w_trig`. This makes
   EVERY existing crossx histogram (pair_pt/eta/minv/dr, pt_150, and the R_AA 3D
   `h3d_op_crossx_..._vs_centr...`) reco-corrected with NO histogram renaming,
   so the crossx plotters and R_AA input require no changes. By construction the
   new nominal == the validated `_corr_unfolded_reco_trig` stage. → `/review-analysis-code`.
3. Rerun crossx RDF (pp + 3 PbPb yrs, skip ntuple).
4. Rerun nominal crossx plotters (`plot_single_b_crossx_{pp,pbpb}.cxx`). → `/review-plot`.
5. R_AA: `RAA_plotting.cxx` is **stale legacy** (Mac `base_dir`, old
   `*_single_b_ana_hists*` filenames, reads an `h3d_ss_crossx_...` the pipeline
   does NOT produce) → not runnable on the cluster; modernization = task_06.
   The reco placeholder is nonetheless now in the R_AA *input* histogram
   (`h3d_op_crossx_..._vs_centr...`). Document; do not silently rewrite R_AA infra.
6. Update placeholder.md / roadmap.

Design note (changes earlier invariant): the §5(d) "trig+reco == existing
trig-only when w_reco=1" invariant no longer applies — by user request the
nominal corrected weight now INCLUDES reco. Trig-only-no-reco nominal is
preserved only in the dated backup; raw/reco/reco+trig stages remain for
before/after visualization.

---
### Progress detail (append-only history)
**Step 3 — DONE (2026-06-15):** `/review-plot` reviewer returned PASS (0
issues) on retry after API recovery. Log:
`.claude/logs/review-plot-20260615-032109-reco-eff-placeholder-digitize.md`.

**Steps 4–5 DONE (code), review next.** Implemented:
- `RDFBasedHistFilling/CorrectionStages.h` (new) — `CorrectionType` /
  `CorrectionStage` enums + `CrossxCorrectionStages()` table.
- `RDFBasedHistFillingData.{h,cxx}` — `OpenRecoEffPlaceholderFile()` +
  `EvaluateSingleMuonRecoEffPlaceholder(centrality,pt,q_eta)` (TGraph store
  `s_reco_eff_ph_map`; pp barrel/endcap, PbPb F.2 ctr mapping; pT clamp [4,19],
  eff clamp (0,1]).
- `RDFBasedHistFillingPP.cxx` FillHistogramsCrossx — reco columns + w_unfold +
  `cw_*` stage columns + 4 staged `h2d_crossx_pair_pt_pair_eta_..._corr_*` hists.
- `RDFBasedHistFillingPbPb.cxx` FillHistogramsCrossx — same, per centrality
  (`h2d_op_crossx_..._vs_pair_eta_vs_pair_pt_<ctr>_corr_*`).
- ACLiC: both `.so` rebuilt clean (only pre-existing ParamsSet sign-compare
  warnings). Standalone key/interp test PASSED (central turn-on 0.55→0.86,
  centrality bump applied, pp barrel 0.92/endcap 0.85).

**Steps 4–5 DONE & REVIEWED (2026-06-15):** `/review-analysis-code` PASS
(iter 1). All 8 specific checks passed; 7 eval anchors MATCH builder. One INFO
(PbPb ctr==-1 → pp graphs, no impact) hardened defensively: PbPb effcy_reco1/2
lambdas return 1.0f when avg_centrality<0. Recompiled clean. Log:
`.claude/logs/review-analysis-code-20260615-234851-reco-eff-placeholder-stages.md`.

**Steps 6–7 DONE (2026-06-15/16):**
- Step 6 backup: plot dirs → `*_backup_20260615` (pp24, pp24_pt_150,
  pbpb_23_24_25_combined[_pt_150]); crossx hist ROOTs →
  `dimuon_data/crossx_hist_backup_20260615/` (pp + 3 PbPb yrs).
- Step 7 rerun: all 4 RDF crossx scripts rc=0 (read existing hadded ntuples).
  Staged hists present: pp 4 (`..._corr_*`), PbPb 6 ctr × {raw,reco,reco_trig}.
  Sanity: bin-by-bin reco ≥ raw in all 581 pp filled bins (0 violations); reco
  inflation ⟨w_reco⟩ ≈ 1.26 (pp) / 1.68 (PbPb ctr10-20). The 30 pp bins with
  reco+trig<reco are the PRE-EXISTING trigger-gap behavior (w_trig=0 on failed
  trig lookup, same as existing crossx_weight_trig_corr), not the reco change.

**Step 8 — DONE (2026-06-16):** new before/after plotter
`plotting_codes/single_b_analysis/plot_crossx_reco_eff_stages.C` → 3 lines
(raw / reco-only / reco+trig) per pair-eta bin, PbPb combined + pp24. Plots in
`dimuon_data/plots/sanity_check_crossx/{PbPb_combined,PP_2024}_reco_eff_stages_pair_pt_in_eta.png`.
`/review-plot` PASS (iter 1; integrals MATCH, physics ordering + pp-vs-PbPb trig
gap confirmed). Log:
`.claude/logs/review-plot-20260616-000444-reco-eff-stages-beforeafter.md`.
**Step 9 (docs) — DONE.** placeholder.md item 3, roadmap ledger updated.

---
> NOTE: the PAUSE/NOT-STARTED block below is **superseded** — all steps 3–9
> completed after API recovery (see COMPLETION SUMMARY at top). Kept for history.

### (prior) Step 3 PAUSE record
**Step 3 — PAUSED (infrastructure, 2026-06-15):** `/review-plot` reviewer
subagent could not be spawned — API returned 529 Overloaded on 5 consecutive
attempts (~15 min). Executor work is complete & self-verified (0–10%
reproduction matches F.2 p01). Resume = re-run `/review-plot` with the same
args (reproduction in `plots/reco_effcy/placeholder/`; originals in
`/tmp/recoeff_pages/`) once the API is healthy.

**Steps 4–9 NOT STARTED.** They are reviewer-gated (`/review-analysis-code`,
`/review-plot`) and equally blocked by the API outage. Resume order: finish
Step 3 review, then Step 4 (correction-stage enum in PP+PbPb RDF), Step 5
(ε₁·ε₂ reco weight + 1/ε clamp), Step 6 (back up current crossx plots), Step
7 (rerun crossx pp+PbPb, SKIP ntuple), Step 8 (3-line before/after plots),
Step 9 (placeholder.md + docs + roadmap ledger).

> Note: `/tmp/recoeff_pages/` is tmp — re-render from PDFs if cleared:
> dimuon note F.2 pp.97–103, HF R_AA note pp.41–42 (see Context paths).
