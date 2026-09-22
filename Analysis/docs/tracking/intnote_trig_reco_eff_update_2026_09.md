# Internal note: trigger-efficiency + reconstruction-efficiency sections — update to the September 2026 procedures (+ figure sync)

**Status:** DONE 2026-09-18 (kept in Active until the user has read §Inconsistencies)
**Started:** 2026-09-17
**Mode:** Implementation (academic writing)
**Predecessor:** `intnote_trig_reco_eff_sections.md` (DONE 2026-08-24 — wrote both sections
from the code of that date). This doc brings them up to the code of 2026-09-17.

## Objective

1. Bring `IntNotes/tex/trigger_efficiency.tex` and `IntNotes/tex/reconstruction_efficiency.tex`
   into agreement with the CURRENT code: the selection change (muon pT 4.5 GeV, gap window
   (−1.25,−1.05), pair pT 9 GeV, pair-level |η^pair| < 2.2, all-vertex pp pairs), the HYBRID pp24
   pair trigger weight (ΔR procedure on the |η^pair| fold below 74.24 GeV; single-value MC pair
   efficiency × SF₁SF₂ above), the rebuilt Step-3 ΔR fits (constrained expo / polynomial,
   baseline screens, interp fallback, no raw-bin tier), the rebuilt pair reco-efficiency map on
   the new selection, and the Pb+Pb changes that touch these sections (signal-region migration,
   centrality acceptance, 2026 data).
2. Re-sync every figure the two sections use (Gate G4), replacing stale ones and adding the
   figures the new procedure needs (single-value pair efficiency, |η|-fold ΔR fits, new stage plot).
3. **Ground truth = the CODE.** Where a tracking doc disagrees with the code, record the
   disagreement in §Inconsistencies for the user; do not silently follow the doc.

## Autonomy Contract (DONE 2026-09-18 — items 1–7 met; item 5 = review loop closed at the cap with zero CRITICAL)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. `tex/trigger_efficiency.tex` describes the delivered pp24 HYBRID weight (§2 of
     `pp24_trig_eff_hybrid_application.md` as verified against the code), the current T&P
     selection (4.5 GeV, gap windows, |η^pair|), the current Step-3 fit family and screens, the
     current closure status, and the current Pb+Pb state; every number re-read from the current
     ROOT/text outputs.
  2. `tex/reconstruction_efficiency.tex` describes the rebuilt pair reco efficiency on the new
     selection (numbers re-read from the 2026-09-10 `pair_reco_eff_pp24_full.root` and its
     value table), the current Pb+Pb placeholder state, and the current signal-region agreement.
  3. `figures/figure_manifest.yaml` updated; `/sync-note-figures` CLEAN (0 stale / 0 missing /
     0 orphaned); new figures added where the procedure changed.
  4. `/compile-note` clean (0 undefined refs/citations, no missing figures).
  5. `/review-note` PASS on both sections.
  6. §Inconsistencies (code vs tracking docs) written and reported to the user.
  7. Committed (IntNotes submodule, then parent), by explicit path.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

Authoritative for this doc: assembled from the CODE. Reference procedures in the parent docs
(`pp24_trig_eff_hybrid_application.md` §2–4, `mc_trigeff_single_value_pair_eff.md` §2–3,
`mc_trigger_efficiency.md` §2–3, `mu_pt45_gap125_pairpt9_adoption.md` §2–3) are cross-checked
against the code in Step 1; any divergence goes to §Inconsistencies.

### 1. Motivation
Every pair in the cross-section carries `w = 1/(ε_trig^pair · ε_reco^pair)`. The note must
describe both factors as they are computed today.

### 2. Top-level equation (pp24, from the code as of 2026-09-17)
```
ε_trig^pair = ε^data_1 ε^data_2 ε_ΔR(ΔR; pT-bin, |η^pair|-group)        9 ≤ pT^pair < 74.24 GeV
            = ε^pair_MC(pT-bin, |η^pair|-group; OS, signal window) · SF_1 · SF_2   74.24 ≤ pT^pair < 150
SF_i = ε^data_i / ε^MC_i ;  ε_trig^pair capped at 1
ε_reco^pair(pT^pair, η^pair, ΔR) from pp24 fullsim, 3D cell → ΔR-integrated → inclusive fallback
```
Pb+Pb: bare mu4 union `ε1+ε2−ε1ε2`, Run-2 single-muon reco placeholder — to be re-verified in Step 1.

### 3. Step-by-step
- 3a. Data T&P single-muon turn-ons (selection now 4.5 GeV, gap windows (−1.25,−1.05),(−0.10,0.06),(2.20,2.40)).
- 3b. MC ε_ΔR on 6 pair-pT bins × 3 |η| groups, OS, expo primary (poly in 3 named forward cells), interp fallback.
- 3c. Single-value MC pair efficiency in bins 7–8 × 3 |η| groups, × SF₁SF₂; D1 pT-merged fallback for one cell.
- 3d. Pair reco efficiency 8 × 9 × 4 on the 9-GeV axes, gap + |η^pair| fiducial on both sides.
- 3e. Pb+Pb: union weight, Run-2 placeholder, signal region now migrated to the pp form.

### 4. Negative constraints
- Do NOT retype binnings: pair-pT from `ParamsSet::pair_pt_coarse_bins` (9→150, 8 log bins), pair-η
  from `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`, gap windows from `ParamsSet`.
- Do NOT present Pb+Pb reco efficiency as a measurement (Gate G7).
- Do NOT edit sibling tracking docs to "fix" a disagreement — record it.

## Implementation Plan

1. **Ground truth gathering (delegated, parallel, read-only; one scratch doc each):**
   - `_sub_note_trigeff_1.md` — trigger chain from CODE + numbers from current outputs + the
     sentence-by-sentence stale list for `trigger_efficiency.tex` + doc-vs-code disagreements.
   - `_sub_note_recoeff_2.md` — same for `reconstruction_efficiency.tex`.
   - `_sub_note_figs_3.md` — figure inventory: current sources' mtimes vs the note copies, new
     figure candidates, manifest edits; other note tex/doc files (datasets.tex, placeholder.md,
     analysis_overview.md, systematic_uncertainties.md) that the sections cross-reference.
2. Rewrite `trigger_efficiency.tex` (main agent). Reviewer: `/review-note`.
3. Rewrite `reconstruction_efficiency.tex` (main agent). Reviewer: `/review-note`.
4. Manifest + `/sync-note-figures`; new plots only if a needed figure does not exist (then `/review-plot`).
5. `/compile-note`; `/review-note`; `/verify-citations` on new cites.
6. §Inconsistencies → user; commit submodule then parent; INDEX.

## Progress Log
(append-only)

- **2026-09-17, Step 0 (done).** Doc triage: read `INDEX.md`; read fully
  `intnote_trig_reco_eff_sections.md`, `pp24_trig_eff_hybrid_application.md`; read the Physics
  Procedure + status sections of `mc_trigeff_single_value_pair_eff.md`, `mc_trigger_efficiency.md`
  (§1–2, R35), `mu_pt45_gap125_pairpt9_adoption.md` (§PP, Remaining Work, Latest Stage). Read both
  current tex sections and the manifest. Verified in `ParamsSet.h`: `pair_pt_coarse_bins` = 8 log
  bins 9→150 (edges 9, 12.79, 18.18, 25.85, 36.74, 52.23, 74.24, 105.53, 150), `pT_bins_150` = 16 log
  bins 9→150, `single_mu_fiducial_gap_cuts = {(−1.25,−1.05),(−0.10,0.06),(2.20,2.40)}`,
  `pair_eta_fiducial_max = 2.2`, `signal_pair_pt_min = 9`, `signal_minv = (1.08, 2.9)`,
  `single_mu_pt_coarse_bins = {4.5,8,14,25,100}`. Noted: INDEX scope of `muon_gap_cuts_acceptance.md`
  still says the barrel/endcap window runs to −1.30 (code: −1.25) → §Inconsistencies candidate.

- **2026-09-18 00:50, user instruction (mid-turn):** after the trigger update, write the reco
  procedure for the gap-acceptance definition as well (plots to follow). Then found that the
  concurrent session had already FINISHED (its Latest Stage: DONE 2026-09-17; new product
  `reco_eff/pair_reco_eff_pp24_full.root` 20:39, applied plots 20:55, pp24 crossx 20:46, R_AA
  21:03) → the reco section was written against the NEW product with its numbers.
- **2026-09-18, Steps 2–5 DONE.**
  * `tex/trigger_efficiency.tex` rewritten (hybrid weight §Strategy + new §sv with the six-cell
    table; T&P 4.0 GeV probe threshold stated; q·η/gap edges 2.2/−1.25; 18 fold cells, shape
    constraints, two baseline screens, poly-primary cells, interp fallback; closure = honest
    PHbox, figure dropped; application 1.88 / seam / sawtooth; PbPb sentinel + signal-region
    PHboxes replaced by the resolved statements; Run-3 trigger citation `TRIG-2022-01` added to
    `ANA-HION-2023-07-INT1.bib` and used).
  * `tex/reconstruction_efficiency.tex` rewritten for the gap-acceptance definition
    (`pair_reco_eff_gap_acceptance.md` §2): truth signal region without gap cuts, reco leg carries
    the data selection; 0.817 acceptance factor; migration 0.02 % in / 3.2 % out-of-range;
    inclusive 0.585 T / 0.636 M (verified from num/den, matches the stored TParameter); vs ΔR
    0.642/0.584/0.513/0.268; vs pT 0.533→0.67; vs η with the two gap dips 0.45/0.48; census
    144 empty + 1 measured zero; reco factor 1.688 on the 2026-09-17 20:46 crossx (trigger factor
    unchanged 1.880, total 3.23); R_AA mixed-fiducial PHbox (0.82 integral, 0.45–0.48 in the gap
    panels, from the concurrent doc's R3 — NOT independently re-derived here).
    Old section backed up to the session scratchpad.
  * Manifest: 16 entries (−closure, −before/after; +singles MC vs data, +ε_ΔR bin-6 cells,
    +reco vs η; ε_ΔR inclusive repointed to `paireta_merged`); two orphans deleted;
    `sync_note_figures.py` → SYNC: CLEAN.
  * Build (CVMFS TeX Live, latexmk): exit 0, 37 pages, 0 undefined refs/cites, 0 overfull,
    0 LaTeX warnings (one "float too large" fixed by 0.98→0.70\textwidth on the 4×3 singles grid).
  * `docs/placeholder.md` rows 3 and 7 + header date updated (item-3 detail had already been
    updated by the concurrent session).

- **2026-09-18 01:15–02:20, Step 6 — `/review-note` loop (5 iterations, log
  `.claude/logs/review-note-20260918-011450-trig-reco-eff-note-update.md`).** Five independent
  reviewers re-extracted every quoted number from the ROOT/text products. Substantive catches,
  all fixed: (1) I had quoted the two dip-row EFFICIENCIES (0.45/0.48) as the gap-ACCEPTANCE factor
  — the factor is new/old = **0.63 / 0.71** (0.81 outermost, 0.96–1.00 elsewhere; recomputed
  myself, `scratchpad/accrow.C`, identical at both WPs); (2) the T&P conditions are muon-level
  (mu4) and pair-level (2mu4) trigger MATCHES within ΔR < 0.02, not event decisions; (3) the
  Step-3 forward low-pT veto (pT < 7 at q·η < −2, MC only) was missing; (4) the kinematic bound
  was inverted — 2m/pT is the MINIMUM opening angle, the maximum is m/√(pT,min(pT−pT,min)) =
  0.64 at 9 GeV, 0.14 at 100 GeV; (5) the closure code is the all-bin cascade incl. a raw tier
  and knows no 74.24 GeV split; (6) poly cells: expo χ²/ndf 4.4/5.7/2.8 (p at bound only in the
  first) vs poly 2.3/2.1/1.5; bin-6 forward baseline error 0.14 % (step-like fit); wide/sig
  −3…+16 %; SF statements qualified (plateau vs turn-on; −2.4<q·η<−2.0 MC fit undershoots its
  saturated points); 12/9/4/5 Pb+Pb fits per year with c railed at 0.1; 6→19 fits > 1 per year;
  the Run 2 placeholder is NOT "measured away from the gaps" (Run 2 had no gap cut) — reworded.
  Loop closed ESCALATED at the 5-iteration cap with ONE trivially-fixed WARNING applied after the
  last review (the placeholder sentence) and not re-reviewed; zero CRITICAL outstanding. Final
  build: 38 pages, 0 undefined refs/cites, 0 overfull; sync CLEAN. `/verify-citations` was run
  inside the loop (iteration 1: 7 keys resolve, 12/13 claims VERIFIED, the MINOR_DISTORTION fixed).

## Inconsistencies (doc vs code) — TO ELEVATE TO USER
(append-only)

Confirmed against the code by the subagents and by me; **no sibling doc was edited.** Code = ground truth.

1. **`muon_gap_cuts_acceptance.md` F17 (+ its INDEX.md scope line)**: barrel/endcap window runs to
   **−1.30**. CODE: `ParamsSet.h:639` = **(−1.25, −1.05)** (changed 2026-09-08). F17's ε_acc
   0.8789/0.8765 and per-window costs are for the superseded window at pT > 4.
2. **`Analysis/docs/placeholder.md` item 7**: "expo / opposite sign / last-two-pair-pT-bins-merged;
   inclusive 0.818 at ΔR→0; `RDFBasedHistFillingPbPb.cxx:1008`". CODE: mode `nocorr_etamerge` (|η^pair|
   fold, no pT merge) + poly cells + interp fallback below 74.24 GeV, single-value × SF₁SF₂ above; the
   inclusive OS curve now starts at 0.780; the Pb+Pb union is at `:1060`. **Item 3**: inclusive
   0.6620/0.7301, ε_acc 0.9133, figure dir `tight/applied/` — all pre-4.5 GeV (current file: 0.7157/0.7780;
   no ε_acc; and the definition is being changed again). **Item 6**: "PbPb 2026 not yet in skim" — the 2026
   lumi/GRL were registered 2026-09-10 and the 2026 crossx/turn-ons exist. (Registry Gate G7 checks the
   note against this file: items 3 and 7 are updated by this task; item 6 is flagged only.)
3. **`mc_trigger_efficiency.md`** §3.0(c)/(c′): "pT > 4 GeV", "truth pT > 4 GeV" — CODE 4.5
   (`MCTrigEffPairSelection.h:101-105,127-129`). Round-8 contract: "8 log bins 8→150" — CODE 9→150.
   R22: forward window {2.30,2.40} applied — CODE {2.20,2.40}. §2 states the Pb+Pb union WITH
   ε_ΔR^single/ε_ΔR^cross as the top-level equation — CODE applies the BARE union.
4. **`mc_trig_eff_closure.md`** Latest Stage: "the delivered correction — including its fallback chain — is
   the object tested". CODE: the crossx uses `DrCorrectionCrossxEvaluator` + `PairTrigEffCrossxEvaluator`
   (hybrid), which `FillMCTrigEffClosure.cxx` does not exercise; all closure numbers are on the retired
   configuration; the raw tier it describes is unreachable from the crossx.
5. **`mc_trigeff_single_value_pair_eff.md`** Latest Stage (2026-09-09): "every number this doc delivers is
   stale … the guard refuses the delivered ROOT file". FILES: `pair_trig_eff_pp24_full.root` was regenerated
   2026-09-10 16:12 on the 9 GeV axis and IS consumed by the crossx. `PairTrigEffEvaluator.h:49-50` and the
   file's provenance still say "NOT wired into the cross-section" — stale since 2026-09-17.
6. **`pp24_crossx_rerun_2026_08.md`** (CLOSED) Physics Procedure gives ε₁ε₂ε_ΔR for the whole spectrum and
   INDEX quotes 0.6620/0.7301/0.9133/×1.477 as key numbers with no STALE marker — superseded by the hybrid
   and the 4.5/9 GeV rerun.
7. **`pp24_all_vertex_pairs.md`** Latest Stage: "MC half NOW RUNNING (2026-09-08)"; "Pb+Pb still has neither
   the fiducial nor the pair-level cut, 44-bin axis". CODE/OUTPUTS: MC half superseded by the 2026-09-10
   rerun; Pb+Pb migrated (`RDFBasedHistFillingPbPb.cxx:977-995`); shared 48-bin axis. Listed ACTIVE.
8. **`systematic_uncertainties.md` §1a** "plateau normalisation … the fit stage reads the plateau
   directly". CODE: the applied `nocorr*` modes never divide by the measured plateau (free baseline C);
   the plateau is only the reference of the `DrCorrBaselineConsistent` screen (|ln(C/plateau)| ≤ ln 1.5).
9. **`muon_wp_registry.md` §4**: placeholder builder writes only Medium, Tight keys hard-coded, pp not on
   `PairRecoEffEvaluator` — CODE: both WPs built, WP-matched keys (`RDFBasedHistFillingData.cxx:813`), pp on
   the 3D pair map. **`reco_eff_placeholder_run2.md`** Latest Stage + INDEX: "pp = HION-2019-58 Fig. 31" —
   pp has used the fullsim 3D map since 2026-08-18.
10. **`analysis_overview.md`**: gap windows and Pb+Pb migration AGREE; only "(Pb+Pb histograms must be
    refilled before R_AA is quoted)" is stale (refilled 2026-09-10 and 2026-09-17).
11. **`mc_trigeff_dr_binning_approaches.md`** 2026-09-08 paragraph "crossx application still
    `nocorr_ptmerge`" survives below its own 2026-09-17 header — CODE `nocorr_etamerge`.
12. **Code comments (not docs)**: `SingleMuEffEvaluator.h:45-46` (q·η contiguous to 2.3 / gap > 2.3 — now
    2.2); `CommonEffcyConfig.h:120-122` + `RDFBasedHistFillingPythiaFullsim.cxx:244-248` ("max truth ΔR
    0.701" — current file < 0.65); `build_pp24_fullsim_pair_reco_eff.C:26-29` (144 empty incl. one
    den-only — now 145 empty, 0 den-only); `CorrectionStages.h:19-22` (w_reco = Run-2 placeholder — pp is
    the 3D map); `ParamsSet.h:313-319` "STILL NOT APPLIED: PbPb crossx, overlay reco-eff, truth
    acceptance, template fits" — all four apply the gap set now.
13. **`.claude/CLAUDE.md` §Binnings** still quotes `pair_pt_coarse_bins` = 8/15/27/50/150 and a
    `_pt150` variant — CODE: 8 log bins 9→150, variant deleted.
14. **INDEX.md** scope of `pbpb2026_analysis_support.md` says step 14 RUNNING; the doc says DONE, awaiting
    experts. (Scope-line staleness; not physics.)
15. **The 2026-08-24 note text itself** carried two Pb+Pb PHboxes (silent-drop sentinel; pp/Pb+Pb signal
    regions differ) that the code has since RESOLVED (fiducial migration 2026-09-08 + centrality
    acceptance filter) — `pp_trig_eff_highpt_jump.md` still reads BLOCKED in INDEX although its item (2)
    is closed by the adoption doc's option (d).
16. **`pair_reco_eff_gap_acceptance.md` R3** quotes "×1/0.45 in |η^pair| ∈ [1,1.5], ×1/0.48 in
    [−0.5,0.5], ×1/0.79 elsewhere" as the per-panel acceptance correction. Those are the NEW
    efficiencies; the acceptance factor (new/old) is 0.63 / 0.71 / 0.81 (outer) / 0.96–1.00
    (`scratchpad/accrow.C`, both WPs). The integral 0.82 is right. (Found by the note reviewer.)
17. **Code comments (found by the reviewers):** `CommonEffcyConfig.h:121` states the opening-angle
    bound as 2m/pT (that is the minimum, not the maximum); `RDFBasedHistFillingData.cxx:780,826`
    say the Pb+Pb placeholder is the "Run 2 Medium fits" (both WPs built, WP-matched key);
    `build_pp24_fullsim_pair_reco_eff.C` header "144 of 288 deliver no efficiency" (145 = 144
    empty + 1 measured zero); `PairTrigEffEvaluator.h:49-50` + the pair-eff file provenance
    "NOT wired into the cross-section" (wired since 2026-09-17).
18. **`IntNotes/tex/datasets.tex` / `introduction.tex`** describe Pb+Pb 2023–2025 only, while
    the 2026 sample is in the combined crossx and has its own turn-ons (`pbpb2026_analysis_support.md`
    defers the chapter update to the note gate chain; two runs under expert review). The trigger
    section carries a \PHtext for this; the datasets chapter is out of this task's scope.

## Results & Observations

### R1 — Trigger chain, ground truth from the code (merged from `_sub_note_trigeff_1.md`, verified 2026-09-17)
- **Equation as coded** (`Utilities/PairTrigEffCrossxEvaluator.h:179-214`, `RDFBasedHistFillingPP.cxx:379-410`):
  region A 9 ≤ pT^pair < 74.24 GeV: ε_trig = ε₁ε₂ ε_ΔR(ΔR; cell), cell = 6 coarse pair-pT bins × 3 |η^pair|
  groups (<1, [1,2), [2,2.2)); ε_ΔR = f/C for ΔR < 1, 1 above; expo primary (A ≤ 0, p ≥ 1, λ ∈ [0.02,3]),
  polynomial primary in pT bins 2–4 × [2,2.2), linear-interpolation fallback, NO raw-bin tier (throw).
  Region B 74.24 ≤ pT^pair < 150: ε_trig = ε^pair_MC(cell; OS, 1.08–2.9 GeV, ≥ 50 raw pairs) × (ε₁/ε_MC,1)(ε₂/ε_MC,2);
  cell [105.5,150)×[2,2.2) (3 raw pairs) served from the pT-merged cell 0.441 ± 0.096 (62 pairs) — D1 TEMP;
  SS pairs use the OS numbers — D2 TEMP. Cap ε_trig ≤ 1 (counted; 0 pairs capped). pT^pair ≥ 150: ε₁ε₂.
  ε_i = data T&P erf×log TF1 at exact pT, cap 1, floor 0.01 (`RDFBasedHistFillingData.cxx:694-746`);
  ε_MC = MC erf×log TF1, pT clamped into [4.5,60], cap 1, floor 0.02 (`SingleMuEffEvaluator.h`).
- **Pb+Pb**: bare union ε₁+ε₂−ε₁ε₂ (`RDFBasedHistFillingPbPb.cxx:1060`), no ε_ΔR; signal region MIGRATED
  2026-09-08 to the pp definition (`:976-995`); centrality-acceptance filter 0 ≤ c < 80 (`:1009-1041`);
  the −1 sentinel/silent-drop path is dead code (coarse q·η binning contiguous over the fiducial region;
  lookup throws instead). A 2026 Pb+Pb turn-on exists (fit file 2026-09-17 05:00).
- **T&P**: tag = HLT_mu4 match, probe = other muon; numerator = event HLT_2mu4; ΔR > 0.8; probe-only gap
  cut; both assignments and both signs summed; per probe charge. T&P NTuple keeps muon pT > 4.0 GeV by
  design (`DimuonDataAlgCoreT.c:628`; nominal is 4.5); data fit range [4,60] with pivot 4; MC fit [4.5,60].
  q·η binning: 10 intervals {−2.4,−2.0,−1.5,−1.0,−0.5,0,0.5,1.0,1.5,2.0,2.2}; gap windows
  [−1.25,−1.05], [−0.10,0.06], [2.20,2.40].
- **pp24 Tight turn-ons** (`single_mu_effcy_pT_fit.root` 2026-09-09 23:38): P 0.683–0.949 (barrel 0.683–0.757),
  m 2.73–5.05 GeV, s 0.62–1.66 GeV, c ≤ 0.073, χ²/ndf 0.52–1.89 (median 1.07; all < 1.5 but one);
  2/20 fits exceed 1 inside [4,60] (both [2.0,2.2), max 1.009). Pb+Pb 2023 Tight: normFermi 0.60–0.952
  (3 railed at 0.6), χ²/ndf 0.36–1.56, 6/120 exceed 1; centrality means normFermi 0.723 (0–5 %) → 0.816
  (50–80 %), ε(60 GeV) 0.816 → 0.881, f(4.5) 0.40 → 0.56. 2024/25/26: normFermi 0.608–0.956.
- **Region-A cells** (OS Tight etamerge fit files 2026-09-10 14:51–14:53): 15 expo + 3 poly + 0 interp;
  delivered ε_ΔR ∈ [0.2837 (bin 6 × [1,2)), 1.000]; barrel ε_ΔR(0) 0.862/0.773/0.626/0.464/0.337/0.286 for
  bins 1–6; expo χ²/ndf 0.62–3.60 (median 1.11); C/plateau ∈ [0.966,1.277]; median |plateau−1| 0.57 %
  (max 17 %), median rel. stat 1.2 %. Poly cells: f(0)/C = 1.000/1.000/0.861, non-monotonic dip to
  0.82/0.77/0.66 near ΔR ≈ 0.3. Inclusive OS expo: C 0.9975, ε(0)=0.780, ε(0.3)=0.959, ε(0.5)=0.999,
  χ²/ndf 126.1/16 = 7.9; inclusive SS: ε(0)=0.307, χ²/ndf 2.07.
- **Single-value cells** (`pair_trig_eff_pp24_full.root` 2026-09-10 16:12; OS, sig, Tight):
  bin 7 [74.2,105.5): 0.222±0.016 [1429/344] / 0.554±0.033 [593/309] / 0.444±0.098 [59/27];
  bin 8 [105.5,150): 0.227±0.034 [294/64] / 0.431±0.085 [101/37] / 0.333±0.272 [3/1] (refused);
  ptmerge [74.2,150)×[2,2.2): 0.441±0.096 [62/28]. Control bin 6: 0.318/0.610/0.607. wide/sig: +13/+10/+12 %
  (bin 7), −3/+15/+64 % (bin 8). SS raw counts bin 7: 33/6/0; bin 8: 6/0/0.
- **SF** ε_data/ε_MC from the current fit files: barrel 0.78–0.87, endcaps 0.96–1.01.
- **Applied on pp24 data** (`histograms_real_pairs_pp_2024_2mu4_nominal.root` 2026-09-17 02:10):
  trigger-only factor 1.880 (mean ε 0.53), reco 1.380, total 2.636; per fine bin trig 2.36 [9,10.7) →
  1.74 [18.2,21.7) → 2.84 [62.3,74.2) → 3.77 [74.2,88.5) → 4.59 [105.5,125.8) → 4.11; per (fine pT × 0.1 η)
  cell 1.27–7.25; region B holds 0.031 % of the raw yield; ⟨w_trig⟩ 1.88 (A) / 3.90 (B).
- **Closure**: NO closure exists on the current inputs or on the hybrid cascade (last fill 2026-09-08
  01:26, pre-4.5 GeV configuration: 0.9951 all-OS / 0.9901 signal, ⟨r₁r₂⟩ 1.279/1.272); the 2026-09-10
  attempt threw on a stale file. Not quotable as a test of the delivered correction.

### R2 — Reconstruction chain (merged from `_sub_note_recoeff_2.md`; DEFERRED per user decision)
- Current product `reco_eff/pair_reco_eff_pp24_full.root` (2026-09-10 14:10, fiducial definition):
  inclusive 0.7157 Tight / 0.7780 Medium; vs ΔR 0.775/0.716/0.634/0.337; 143 measured / 145 empty cells;
  max generated ΔR < 0.65; applied reco factor 1.380 on the 2026-09-17 crossx. Axes 9→150 / |η| < 2.2.
  ε_acc has NO current value. Pb+Pb: Run-2 single-muon placeholder product unchanged (clamp [4,19]),
  signal region identical to pp, 2026 on the 2023 FCal calibration; overlay RDF class uses the full gap
  set; only test overlay samples. **All of this is about to be superseded** by
  `pair_reco_eff_gap_acceptance.md` (gap cuts off the truth denominator) — the section is written after
  that product lands.

### R3 — Figures (merged from `_sub_note_figs_3.md`)
- Manifest: 8 STALE by content (sources 2026-09-09…17), 5 MISSING-SOURCE (4 reco applied PNGs; closure),
  2 current (Pb+Pb placeholder). Applied reco PNGs regenerated 2026-09-17 20:27 by me for both WPs
  (then moved by the concurrent session to `before_gap_acceptance/`); will be regenerated again on the
  new product. ε_ΔR figure repointed to `step3_dr_fit/no_plateau_correction_paireta_merged/expo/sign_sepr/`.
  New trigger figures adopted: `step1_singles_data_mc/step1_eff_pt_in_q_eta_bins_muplus.png` (MC vs data
  singles, the SF), `..._pairpt_52.2_74.2.png` (last region-A bin, per-cell fits). Closure figure dropped.

## Remaining Work
- Commit (IntNotes submodule, then parent) — Step 6.
- After the user's decisions: closure of the hybrid on the current inputs (closure code must be
  brought to the hybrid form first); datasets chapter update for 2026; `pair_reco_eff_gap_acceptance.md`
  R3 wording (sibling doc, not edited here).

## Latest Stage
DONE 2026-09-18. IntNotes committed (6207231, 017d048, a8dbb6b); parent commit with this doc, INDEX, placeholder.md and the submodule pointer. Scratch docs merged (R1–R3, §Inconsistencies) and deleted. Open for the user: §Inconsistencies (18 items) and the Remaining Work list.

(superseded) Step 2: rewriting `IntNotes/tex/trigger_efficiency.tex` from R1 (trigger only; reco deferred per the user decision). Figures for it: turn-ons ×2, map, Pb+Pb23 0–5 %, step-1 singles MC vs data, ε_ΔR inclusive (paireta_merged), per-cell bin-6 fits; stage figure referenced from the reco section. Closure figure dropped.

- **2026-09-17, Step 1 (in progress).** Launched 3 read-only subagents: `_sub_note_trigeff_1.md`,
  `_sub_note_recoeff_2.md`, `_sub_note_figs_3.md`. Main agent meanwhile verified
  `Utilities/PairTrigEffCrossxEvaluator.h` (the hybrid, region split = coarse edges N−2, SF per leg,
  cap at 1 counted, D1 single-cell fallback with a throw on a second refusal) and found the Run-3
  trigger reference for the old "citation needed" placeholder: JINST 19 (2024) P06029,
  arXiv:2401.06630, CERN-EP-2023-299 (KB `physics/detector/atlas_run3_muon_performance.md`; PDF
  first page checked) — not in `bib/ATLAS.bib`, to be added to `ANA-HION-2023-07-INT1.bib`.

- **2026-09-17 late, USER DECISION (AskUserQuestion).** A concurrent session
  (`pair_reco_eff_gap_acceptance.md`, opened 2026-09-17) is changing the ε_reco DEFINITION (gap cuts
  off the truth denominator; ε_reco absorbs the gap acceptance; no ε_acc), followed by a fullsim
  refill → builder → pp24 crossx rerun. It moved the applied plots I had regenerated at 20:27 from
  `tight/applied/` to `before_gap_acceptance/{tight,medium}/applied/` (its D3 backup). **User chose:
  "Trigger now, reco after"** — rewrite `trigger_efficiency.tex` now against the current code; write
  `reconstruction_efficiency.tex` only after that session's new `pair_reco_eff_pp24_full.root` and
  crossx refill land (poll that doc's Latest Stage + the product mtime). Consequence for the trigger
  section: quote only reco-independent numbers (trigger-only factor, ⟨w_trig⟩ table, census); the
  before/after and stage figures will be re-synced after the concurrent crossx refill.
- **2026-09-17, Step 1 DONE.** All three subagents returned; scratch docs `_sub_note_trigeff_1.md`
  (300 lines), `_sub_note_recoeff_2.md` (248 lines), `_sub_note_figs_3.md` (149 lines). Merged
  below in §Results and §Inconsistencies.
