# Remove ΔR>0.05 from the single-b signal selection (interim nominal)

**Mode:** Implementation. **Created:** 2026-06-22.
**Reviewer rules:** code → `/review-analysis-code`; plots → `/review-plot`.
**Dependency map (authoritative for the rerun blast radius):**
`Analysis/docs/signal_selection_change_impact.md`.
**Siblings:** `analysis_overview.md` §2 (signal region), `raa_from_rdf_crossx.md`,
`placeholder.md`.

---

## Objective
Remove the `dr > 0.05` cut from the single-b dimuon **signal selection** (data +
truth + fullsim defs), make the no-ΔR-cut spectra the **interim nominal**, and
rerun the full downstream chain (data crossx pp + pbpb 23/24/25, pythia truth
acceptance, all crossx/R_AA/acceptance/cutflow plots). Keep the change reversible
(backup nominal first).

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation
The signal is a single-b dimuon: both muons from one b chain → intrinsically
**small opening angle**, with ΔR ~ 2m/p_T (collimation). The `dr > 0.05` cut was
added only as a **statistics workaround** for the *data-based* dR-dependent
trigger-efficiency inverse-weighting (the `minv>1.08` cut already drives the ΔR
distribution toward zero at small ΔR, starving the small-ΔR bins). It is **not a
physics requirement** of the signal. As shown by the cutflow and the ΔR-vs-pair-p_T
diagnostics (`plots/single_b_analysis/dr_vs_pair_pt_diagnostic/`), a fixed ΔR>0.05
removes a **growing fraction of genuine signal at high pair p_T** (band crosses
below 0.05 by ~40–50 GeV), biasing dσ/dp_T and R_AA(p_T) exactly where b energy
loss is measured.

### 2. Decision
The dR trigger-correlation correction is being moved to **MC** (unbiased trigger
sim). Until that MC exists we cannot determine the correct, possibly
pair-p_T–dependent ΔR cut — and we may need **none**. Interim nominal = **no ΔR
cut**. A ΔR cut may be re-added later (same blast radius, reverse of this doc).

### 3. New signal region (data, reco; pp & PbPb identical)
```
minv > 1.08 && minv < 2.9 && pair_pt > 8
  && m1.charge*m1.eta < 2.2 && m2.charge*m2.eta < 2.2          (NO dr cut)
```
Truth analog: `truth_*` vars + `from_same_b`, no `truth_dr` cut. The ntuple-level
`mindR=0.02` (DimuonDataAlgCoreT) still bounds data from below — this is a
trigger-matching floor, not a signal cut, and is unchanged.

### 4. Negative constraints
- Do NOT touch ntuple processing / `muon_pairs_*` trees, nor the trigger-eff ε^nc
  fits (independent of pair cuts) — see dependency map §1.
- The ΔR-binned diagnostic 2D-hist **axis floors must move 0.05→0.0** so the
  newly-admitted 0.02–0.05 pairs are not lost to underflow; the main crossx
  p_T/η spectra have no ΔR axis and are correct regardless.
- Reco-eff is a Run-2 **placeholder** (pT·q·η, ΔR-independent); editing the
  fullsim `pass_signal_*` defs is code-consistency only — numerically inert for
  nominal until real Run 3 MC. Do the edits; rerunning fullsim is optional now.

## Implementation Plan
1. [in progress] Tracking doc + backup nominal outputs. Register in CLAUDE.md.
2. Code edits via `/review-analysis-code` (quote Physics Proc §3/§4): remove dr
   from the 8 active files + cutflow `kCuts`; extend ΔR axis floors 0.05→0.0
   (PP 457/519/534; PbPb 1000/1073/1089/1136/1154). Files per dependency map §2.
3. Recompile (ACLiC); rerun data crossx pp24 + pbpb23/24/25; rerun pythia truth
   (3 modes). (Fullsim/overlay reco-eff: edit only, defer rerun — placeholder.)
4. Replot via `/review-plot`: crossx pp + pbpb combined, sanity/stage, R_AA
   (mode 6), signal acceptance + cutflow, dr-vs-pt diagnostic refresh.
5. Update docs: analysis_overview §2 (drop ΔR>0.05 from signal region, note
   interim), roadmap, placeholder.md, dependency map (mark executed), status.
6. Verify: crossx/R_AA shapes — expect high-pT yield/acceptance to RISE (recovered
   small-ΔR signal); acceptance cutflow ΔR step should vanish.

## Progress Log
- 2026-06-22 — Step 1: doc created; nominal outputs backed up
  (`dimuon_data/crossx_hist_backup_20260622_pre_dr_cut_removal/` for pp+pbpb crossx;
  `pythia_truth_full_sample/hist_backup_20260622_pre_dr_cut_removal/` for truth);
  user approved full removal + rerun. Registered impact doc as must-read in
  CLAUDE.md Documentation References, analysis_overview.md §2, README.md.
- 2026-06-22 — Step 2 DONE (`/review-analysis-code` PASS iter 2, log
  review-analysis-code-20260622-010221-remove-dr-signal-cut.md). Removed
  `dr>0.05`/`truth_dr>0.05` from all active signal-region/pass_signal/template/
  cutflow defs: PP 382/552, PbPb 920/1194, PythiaTruth 412 (acceptance) +
  325-327 (`kin_cuts` template region — MISS caught by reviewer in iter 1, fixed),
  PowhegTruth 158, PythiaFullsim 129/131, PythiaFullsimOverlay 69/71,
  PowhegFullsim 185/186, cutflow kCuts. ΔR diagnostic 2D-axis floors 0.05→0.0
  (PP 457/519/534; PbPb 1000/1073/1089/1136/1154; 50 bins/upper 1.0 kept). All 7
  RDF classes ACLiC-clean. Only retired legacy files (FillHistogramsCrossx_PP_clean,
  SingleBAnalysisBase) still carry the old cut (intentional).
- 2026-06-22 — Step 3 LAUNCHED (background): crossx chain pp24→pbpb23→24→25
  (`/tmp/dr_rerun_logs/crossx_chain.log`) + pythia truth nonprivate 5.36 nocuts
  (`/tmp/dr_rerun_logs/pythia_truth_chain.log`). Fullsim/overlay reco-eff NOT rerun
  (placeholder-inert per §4; code edited for consistency only).

## Results & Observations
- Reviewer iter-1 catch: the truth template-fit region `kin_cuts`
  (`FillHistogramsTemplateMinvSignalRegion`, PythiaTruth ~325) is a SECOND truth
  signal-region def beyond the acceptance `signal_cuts`; both needed the dr removal.
  Recorded for the dependency map (future cut changes must hit both).
- **Step 3 reruns DONE (2026-06-22), rc=0, no errors:** crossx pp24 + pbpb23/24/25
  + pythia truth 5.36. Self-consistency verified vs current input trees: pp OS
  signal entries 788579 = tree no-dr count; pbpb23 OS 3D entries 83800 ≈ tree no-dr
  83818. **New runs are correct.**
- **dr-removal physics effect (data) — SMALL inclusively, SHAPE at high pT:** the
  ntuple mindR=0.02 + steeply-falling spectrum mean removing dr>0.05 adds only
  +172 pp pairs (+0.02%) / +52 pbpb23 pairs (+0.1%) inclusively, but the recovered
  pairs concentrate at high pair pT: pp crossx ×1.13 (pT>60), ×1.27 (pT>80). The
  TRUTH acceptance is the clean demonstration: α flattens to ~0.72 across all pT
  (was 0.57/0.48/0.38 at pT>60/80/100 — the ΔR-collimation loss is gone). Matches
  the cutflow + dr-vs-pt diagnostics. Physically expected; not a regression.
- **FINDING (separate, pre-existing — surfaced to user):** the pp "inclusive yield
  drop vs backup" (5732→5673, −1%) is NOT from the dr cut. The pre-removal backup
  pp nominal was built from the older/larger `muon_pairs_pp_2024_2mu4.root` (3.37M
  OS, 2025-04-17), while the rerun correctly uses the canonical `_mindR_0_02.root`
  skim (3.30M OS, mtime 2026-06-10) — the proper mindR=0.02 trigger-matched input
  (class candidate order picks it first). So the new pp nominal switched to the
  correct input AND dropped the dr cut. The 06-10 pp nominal (and 06-19 R_AA) may
  have used the older plain file — worth a separate check, but the NEW state is the
  more-correct one. PbPb used `single_mu4_mindR_0_02_res_cut_v2.root` (correct).

## Remaining Work
- None for this rework. **Open follow-up (separate):** the pp "inclusive drop"
  finding — confirm whether the 06-10 pp nominal used the older plain
  `muon_pairs_pp_2024_2mu4.root` instead of the canonical `_mindR_0_02.root`
  (the rerun now uses the correct mindR_0_02 input). Worth a quick check next
  session; does not block this rework.
- Changes uncommitted — awaiting user go-ahead to commit.

## COMPLETION SUMMARY (2026-06-22) — doc CLOSED
ΔR>0.05 removed from the single-b signal selection (interim nominal). Steps 2–4
(`/review-analysis-code` PASS iter 2, reruns, `/review-plot` PASS iter 1) +
Step 5 docs done. Validation: pythia truth signal acceptance now FLAT ~0.72 across
all pair pT (was 0.57/0.48/0.38 at pT>60/80/100 — ΔR-collimation loss eliminated);
cutflow has no ΔR step and all-pT/pT>60 curves track; data crossx inclusive
≈unchanged (+0.02% pp, +0.1% pbpb), high-pair-pT yield +13–27%; R_AA physical &
centrality-ordered; crossx ΔR-2D axis now spans 0–1 with the low-ΔR band populated.
Reco-eff/det-resp fullsim defs edited for consistency (placeholder-inert until real
MC). Reusable blast-radius map: `docs/signal_selection_change_impact.md`. Logs:
review-analysis-code-20260622-010221-remove-dr-signal-cut.md,
review-plot-20260622-013539-dr-cut-removal-replots.md.

## Latest Stage
*(cleared — rework COMPLETE; both review loops PASS. Changes uncommitted, awaiting
user go-ahead to commit.)*
