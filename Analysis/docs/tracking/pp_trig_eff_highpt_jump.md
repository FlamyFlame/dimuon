# Investigation: pp24 trigger-efficiency high-pT jump in crossx / R_AA

**Mode:** Investigation. **Created:** 2026-06-20. **Status:** ROOT CAUSE FOUND
(reviewed `/review-investigation` PASS iter 1, log
`.claude/logs/review-investigation-20260620-050718-pp-trig-eff-highpt-jump.md`).
Fix not yet applied — awaiting user decision on approach.

## Objective
Explain the large discontinuity / sudden jump UP in the pp24
trigger-efficiency-corrected cross-section (and R_AA, which drops/spikes at the
same pair pT) around pair pT ~50–60 GeV. Visible in
`plots/sanity_check_crossx/PP_2024_reco_eff_stages_pair_pt_in_eta.png` (red
"Reco×trig" curve shoots ~1–2 orders of magnitude above raw/reco-only at high pT,
huge error bars). No PbPb centrality shows it.

## Context
- pp crossx weight: `RDFBasedHistFillingPP.cxx:400-462` — 2mu4 PRODUCT weight
  `effcy_pair = ε₁·ε₂`, `w_trig = 1/effcy_pair`, **no floor/cap**.
- PbPb crossx weight: `RDFBasedHistFillingPbPb.cxx:933-985` — single-mu4 UNION
  weight `effcy_pair = ε₁+ε₂−ε₁ε₂`, `w_trig = 1/effcy_pair`.
- Single-muon ε lookup: `RDFBasedHistFillingData.cxx:596-625`
  (`EvaluateSingleMuonEffcyPtFitted`): `val = TF1::Eval(pt); if(val<0.01)val=0.01;`
- Fit files: `pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit.root`
  (pp); `pbpb_20YY/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root`
  (PbPb). Fits made separately (P2 trig-eff pipeline).

## Accumulated Findings (append-only)
1. **The turn-on TF1s return EXACTLY 0 above their fit range [4,60] GeV on
   read-back.** They are compiled-function TF1s (no analytic TFormula —
   `GetExpFormula()`==""). When such a TF1 is written to a ROOT file, only the
   sampled grid over [xmin,xmax] persists; `TF1::Eval(pt)` returns 0 for pt
   outside [4,60]. Verified pp TF1 `..._1_60_TO_2_00_divided`: Eval(60)=0.9893,
   Eval(62.5)=0.0000, Eval(70)=Eval(100)=0.000.
2. **Lookup floors ε to 0.01** (`RDFBasedHistFillingData.cxx:606`): the returned 0
   becomes 0.01. So every muon with pt>60 GeV gets ε=0.01 instead of its true
   ~0.97–0.99 plateau.
3. **pp 2mu4 product weight amplifies it, unfloored.** `w_trig = 1/(ε₁·ε₂)`: one
   >60 GeV muon → w_trig≈×100; both → ×10⁴. (Contrast w_reco, which IS floored at
   0.05 → ≤×20.) High pair_pt has tiny raw yield (0.015–1/bin) so a few
   mis-weighted pairs dominate → the jump.
4. **Data confirmation (pp).** trig×reco correction factor (corr/raw) per pair_pt
   bin from `histograms_real_pairs_pp_2024_2mu4_nominal.root`
   (`h2d_crossx_..._w_signal_cuts` / `..._no_trig_corr`, ProjectionX): smooth
   6.14→1.49 up to [48.7,58.3]; then **12.96 @[58.3,69.8], 73.7 @[69.8,83.6],
   80.8 @[83.6,100.2], 88.1 @[100.2,120]**. Onset bin straddles the 60 GeV
   single-muon fit edge.
5. **PbPb is immune via the union weight, NOT better fits.** PbPb fits ALSO return
   0 beyond [4,60] (Eval(70)=Eval(100)=0). But the union `ε₁+ε₂−ε₁ε₂` is bounded
   below by max(ε₁,ε₂): a floored leading muon (0.01) + an in-range subleading
   (~0.98) → denominator ≈0.98 → w_trig≈1.02. PbPb 2023 trig-corrected op yield
   falls smoothly 5762→…→1.86 with NO jump.

## Ruled Out (append-only)
- **Wrong file path / histogram name when applying the correction** — paths and
  hist names are correct (verified the code opens the right fit + hist files and
  reads the right histograms).
- **Obsolete/wrong 2D gap-region efficiency fallback** — not the cause; the gap
  2D fallback (`RDFBasedHistFillingData.cxx:612-619`) is a separate branch and the
  jump is driven by the fitted-TF1 branch returning 0 out of range.
- **Fit failure / wild fitted function** — fits are fine WITHIN [4,60]; the issue
  is read-back extrapolation returning 0, not a bad fit.
- **Breakdown of the no-correlation assumption (hypothesis 3)** — NOT the cause.
  The artifact is a sharp ×9 step localized exactly at the 60 GeV fit edge, not a
  smooth dR-correlation bias. mu4_mu4noL1 also requires both muons (product-type
  weight) using the SAME out-of-range fits → the jump would persist. The proposed
  full-chain mu4_mu4noL1 rerun is therefore unnecessary to explain THIS artifact
  (it can still be done later to study the correlation assumption itself, but it
  would not remove this jump).

## Recommended fix (follow-up — `/review-analysis-code`, NOT this investigation)
Make the lookup not propagate the out-of-range 0 as ε=0.01. Preferred (a):
**clamp pt into the TF1 range before Eval** (`pt → min(pt, xmax)`), so high-pT
muons get the plateau ε (~0.98) — physically correct since the trigger turn-on is
flat above threshold; fixes pp and PbPb consistently. Alternatives: (b) re-save
fits as analytic TFormula TF1s so Eval extrapolates the plateau; (c) floor/cap
w_trig (symptom only). After the fix: recompile RDF, refill pp (and PbPb) crossx,
replot sanity_check_crossx + crossx + R_AA, confirm the jump is gone.

## Fix applied & verified (2026-06-20) — DONE
7. **Fix (approach a) implemented and reviewed.**
   `RDFBasedHistFillingData.cxx::EvaluateSingleMuonEffcyPtFitted` (fitted-TF1
   branch): added `GetRange(xmin,xmax); pt_eval = clamp(pt, xmin, xmax);
   val = Eval(pt_eval)` before the existing `[0.01,1.0]` guard. So muons above the
   [4,60] turn-on range now get the plateau ε (~0.98) instead of 0→0.01. 2D gap
   fallback and the −1.0 return are unchanged; pp/PbPb weight formulas untouched;
   no cap added to w_trig (cause fixed, not symptom). Shared lookup → fixes pp and
   PbPb consistently. `/review-analysis-code` PASS iter 1 (log
   review-analysis-code-20260620-053030-trig-eff-plateau-clamp.md).
8. **Refilled** pp24 crossx + PbPb 23/24/25 (recompiled RDF, rc=0 each).
9. **Verification (corr/raw ratio per pair_pt bin, BEFORE → AFTER):**
   [58.3,69.8] 12.96→1.32; [69.8,83.6] 73.74→1.32; [83.6,100.2] 80.84→1.13;
   [100.2,120] 88.14→1.62. Low/mid-pt unchanged. Jump ELIMINATED. (Reviewer
   independently re-extracted: all MATCH.)
10. **Replotted & reviewed** (`/review-plot` PASS iter 1, log
    review-plot-20260620-054254-trig-clamp-crossx-raa-replot.md): PP stage plot
    red now tracks blue/black smoothly (no jump); pp24 crossx smooth falling
    spectrum; PbPb stage/crossx unchanged; R_AA vs pT/η/ctr smooth, no
    drop-to-zero discontinuity, PLACEHOLDER label present.

## Latest Stage
*(see REOPENED 2026-06-21 below — interim plateau-clamp being replaced by the
proper continuous-function fix.)*

---

## REOPENED 2026-06-21 — proper root-cause fix: persist & use the CONTINUOUS fit function

**Trigger (user):** the 2026-06-20 plateau-clamp was an interim patch — it still
relies on the broken sampled-grid read-back and evaluates at xmax for pt>60 instead
of the true function. Two deeper points:
1. **Root cause confirmed:** the turn-on TF1s are built from a C++ **lambda**
   (`SingleMuEffcyPtTurnOnFitter.cxx:182,225,268`), not a TFormula string. ROOT
   can't persist a function pointer, so `Write` saves only the `fSave` sampled grid
   over [4,60]; read-back `Eval` interpolates that grid and returns 0 outside it.
   This is ROOT behavior triggered by the lambda choice — NOT a runtime decision.
2. **The fix must use the continuous analytic function at the EXACT pt** (never a
   sampled/nearest point), for BOTH pp and PbPb (the "returns 0 above 60" is
   unacceptable for PbPb too even though the union weight hid it). Cap at usage with
   `min(f,1)` (the multiplicative log term can exceed 1 at large pt). General lesson
   recorded in memory `feedback_use_continuous_fit_function`.

### Exact fit forms (from the fitter; pT_min=4, pT_max=60)
- erf+log (pp):    `[2]*0.5*(1+Erf((x-[0])/(sqrt2*[1])))*(1+[3]*log(1+(x-4)/4))`
  pars [0]mean [1]sigma [2]plateau [3]corrCoef.
- fermi+log (PbPb): `[0]/(1+exp(([1]-x)/[2]))*(1+[3]*log(1+(x-4)/4))`
  pars [0]normFermi [1]pT0 [2]Delta [3]corrCoef.

### Plan
1. **Fitter** `SingleMuEffcyPtTurnOnFitter.cxx`: build the TF1s from **TFormula
   strings** (all 4 modes) so the analytic formula persists on Write; keep identical
   par names/init/limits/fit options. → /review-analysis-code.
2. **Re-fit** (pipeline Stage 6, `single_muon_trig_effcy_pT_fitting()`) for pp +
   PbPb 23/24/25 (reads existing fine-q·η eff histos). Verify: GetExpFormula()≠"";
   Eval(70)/Eval(100) give the plateau (not 0); fitted params ≈ old.
3. **Lookup** `RDFBasedHistFillingData::EvaluateSingleMuonEffcyPtFitted`: REMOVE the
   plateau-clamp; `Eval(exact pt)`; cap `min(val,1)`; keep the 0.01 defensive floor.
   → same /review-analysis-code loop.
4. **Rerun crossx (pp + PbPb 23/24/25) + R_AA + sanity/crossx plots** via a SEPARATE
   agent; verify smooth (no jump) and corr/raw ≈ post-clamp values (~1.3 at high pt).
5. **Chain-wide audit** (separate agent): find any other place a fit is resampled
   before evaluation, or a TGraph/interp is read at the nearest sample point instead
   of interpolated. Report; fix any found.
6. Docs + memory finalize.

### Progress (2026-06-21) — COMPLETE
- **Step 1 (fitter) + Step 3 (lookup):** DONE, `/review-analysis-code` PASS iter 1
  (log review-analysis-code-20260621-204824-trig-fit-tformula-persist.md). All 3
  fit constructors rebuilt from TFormula strings; lookup evaluates exact-pt Eval +
  min(.,1) cap, plateau-clamp removed; 2D fallback + −1.0 return untouched. RDF PP
  ACLiC clean; fitter runs interpreted (pre-existing missing-includes, unchanged).
- **Step 2 (re-fit):** pp re-fit verified — 40/40 TF1s persist formula; params
  IDENTICAL to old (4.349/1.005/0.8977/0.03767); in-range Eval matches old exactly
  (60→0.9893); above range now continuous (70→0.9945, 100→1.0066) vs old 0.0000.
  PbPb 23/24/25 re-fit (separate agent): 240/240 TF1s each persist formula;
  Eval(70/100)≈0.96–1.06 (NOT 0); file size 354KB→195KB (TFormula vs lambda grid).
- **Step 4 (rerun):** crossx refilled pp + PbPb 23/24/25 (rc=0); R_AA + crossx +
  sanity stage plots regenerated. **pp corr/raw ratio smooth** — high-pt bins
  1.32/1.32/1.13/1.61 (was 12.96/73.74/80.84/88.14); low/mid unchanged. PbPb23 OS
  yield smooth/monotonic. PP stage plot + R_AA-vs-pT visually smooth (no jump),
  matching prior reviewed-PASS state.
- **Step 5 (chain audit):** DONE — NO other instance of the bug class. PbPb reco-eff
  placeholder = TFormula logistic (exact Eval, safe); pp reco-eff = TGraph::Eval
  (linear interp, correct); 2D gap fallback = GetBinContent(FindBin) on the RAW
  MEASURED hist (not a resampled fit, correct); event-sel cut graphs = TGraph::Eval
  clamped (correct); centrality = lower_bound percentile (discrete by design); ZDC
  EMG lambda-TF1s = in-memory only, never persisted. One defensive note:
  EvaluateSingleMuonRecoEffPlaceholder Evals without a range guard — safe now
  (analytic fits + pt clamped to fit range); add a guard only if those are ever
  rebuilt from compiled functions.
- **Memory:** `feedback_use_continuous_fit_function` (general lesson);
  `project_tf1_eval_out_of_range` updated (proper fix supersedes the clamp).

### COMPLETION SUMMARY (2026-06-21) — doc CLOSED again
Root cause was the fitter building turn-on TF1s from C++ lambdas (ROOT can't persist
a function pointer → Write saves only the fSave grid over [4,60] → read-back Eval
returns 0 outside). Fixed at the source: fits rebuilt as TFormula strings (formula
persists + extrapolates); lookup evaluates the continuous function at exact pt with
min(.,1) cap. pp + PbPb re-fit and full crossx/R_AA chain rerun; high-pT jump gone;
in-range numbers unchanged. Chain-wide audit found no other instance. The interim
plateau-clamp (2026-06-20) is superseded by this exact-pt fix. `/review-analysis-code`
PASS. Uncommitted — awaiting user go-ahead to commit.

## Latest Stage
*(cleared — continuous-function root-cause fix COMPLETE; chain rerun verified;
chain-wide audit clean. Uncommitted.)*
