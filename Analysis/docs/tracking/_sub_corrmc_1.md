# _sub_corrmc_1 — Corrected-MC study (SF = ε_data/ε_MC) for the MC trigger efficiency

Subagent scratch doc (source of truth). Task: round-7 Autonomy Contract item 5 of
`mc_trigger_efficiency.md` — "Corrected-MC study".

## Mandate (from the orchestrator, verbatim intent)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress.
- Done =
  1. Corrected Step-1 fill (numerator weighted by `SF = ε_data/ε_MC(pT, q·η)`), fitted with the
     SAME fitter → `single_mu_effcy_pT_fit_mc_corrected<wp>.root`; plots in `step1_corrected_mc/`
     (+ `medium/`), corrected-MC vs DATA, with the corrected fit curve; **quantified** agreement.
  2. Corrected Step-3 / Step-4 fills (numerator weight `w·SF/ε_corr`), plots as **1 PNG per
     pair-pT bin, 1 subplot per pair-η bin**, original-MC vs corrected-MC overlaid, in
     `step3_corrected_mc/` and `step4_corrected_mc/` (+ `medium/`); round-7 conditional error
     bars carried through; max/typical |corrected − original| per cell in units of the error.
  3. Driver script `Analysis/pipelines/run_mc_trigeff_corrected.sh`.
  4. Samples: pp_full + overlay, both WPs (Tight nominal, Medium).
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment; when unsure
  whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## ENVIRONMENT BLOCKER (found 2026-08-03, first action) — READ THIS FIRST
The git worktree handed to me
`/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/.claude/worktrees/agent-a3b219f0e24ce140d`
is checked out at commit **b3052ca**, which **predates the whole MC-trigger-efficiency program**:
`Analysis/RDFBasedHistFilling/FillMCTrigEffHists.cxx`, `FitMCSinglesEffcy.cxx`,
`Analysis/plotting_codes/trig_effcy/mc_based/`, `Analysis/docs/tracking/INDEX.md` and
`mc_trigger_efficiency.md` **do not exist there**. (Repo master = 3be0c5a, the round-7 commit.)

Action taken (no git run — HARD RULE respected): I copied the CURRENT master content of the
files/dirs I need from the main working tree into the worktree, and do all editing/compiling
there:
- `Analysis/MuonObjectsParamsAndHelpers/` (whole dir)
- `Analysis/Utilities/` (whole dir)
- `Analysis/RDFBasedHistFilling/{FillMCTrigEffHists.cxx, FitMCSinglesEffcy.cxx, CommonEffcyConfig.h}`

**⇒ The orchestrator MUST NOT `git merge` this worktree branch.** Everything outside those
paths in the worktree is a stale b3052ca snapshot and a merge would revert master. Cherry-pick
the file list in the final summary instead.

## Physics design (follows mc_trigger_efficiency.md §3.1/§3.3/§3.4 + §4)

`SF(pT, q·η) = ε_data(pT, q·η) / ε_MC(pT, q·η)`, both from the fitted TF1s evaluated
CONTINUOUSLY at the exact (pT, q·η) (never resample-to-nearest).
- ε_MC: `f_mc_pt_vs_q_eta_<muplus|muminus>_<qeta>` in `single_mu_effcy_pT_fit_mc<wp>.root`
  (per-sample dir) — the existing `MCEffEvaluator`.
- ε_data: `f_pt2nd_vs_q_eta2nd<ctr>_<sign1|sign2>_2mu4_sepr_py_<qeta>_divided` in
  `.../pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit<wp>.root`
  (pp / noovl) or `.../pbpb_2023/trg_effcy_pT_fitting_to_fermi_plus_log/...` with
  `ctr = _ctr0_5` (overlay, doc D2). sign1 = μ⁺, sign2 = μ⁻.
  2D gap fallback: `h_pt2nd_vs_q_eta2nd<ctr>_<sign>_2mu4_sepr_divided` (x = q·η, y = pT) in the
  data `histograms_real_pairs_..._single_mu4_fine_q_eta_bin<wp>.root`.
  Same clamp-to-fit-range / cap 1.0 / floor 0.02 guards as `MCEffEvaluator` (the
  `project_tf1_eval_out_of_range` lesson: a read-back compiled TF1 returns 0 outside [4,60]).

**"Corrected MC" = every MC muon that FIRED carries the extra weight SF; the denominator is
untouched** ⇒ `ε_corr = Σ_fired w·SF / Σ_all w = ε_MC·⟨SF⟩ ≈ ε_data`.

Step 3 corrected numerator weight: `w · SF₁SF₂ / (ε_corr,1 ε_corr,2)`.
Step 4 corrected numerator weight: `w · SF_leg / ε_corr,leg`.
`ε_corr` from the CORRECTED fit file (self-consistency). If ε_corr were exactly ε_data then
`SF/ε_corr = 1/ε_MC` = the ORIGINAL weight ⇒ the ΔR corrections are identical bin-by-bin. That
is the validity test.

**Conditional (round-7) error terms, generalised.** `N = Σ t_i a_i`, `Var(N) = Σ a_i² p_i(1−p_i)`
(+ leg-leg covariance for Step 4), with a_i the NUMERATOR weight and p_i the TRUE Bernoulli
probability, which is governed by ε_MC (SF is an analysis weight, not a change of the trigger):
- `A = Σ_fired a_i²`, `B = Σ_fired a_i² · ε_MC,i` ⇒ `Var = A − R·B`
- Step 4 covariance: `covP = Σ_{both fired} 2 a₁a₂`,
  `covQ = Σ_all 2 a₁a₂ ε_MC,1 ε_MC,2` ⇒ `+ covP − R²·covQ`
Setting `a_i = w/ε_MC,i` (SF ≡ 1, ε_corr ≡ ε_MC) reproduces the existing nominal expressions
exactly (`A = Σw²/ε²`, `B = Σw²/ε`, `covP = 2Σw²/(ε₁ε₂)`, `covQ = 2Σw²`). Verified algebraically.

## Progress log

### Step 0 (done) — background reading + environment
- Read `mc_trigger_efficiency.md` §1, §2, §3.0–§3.4, §4, round-7 Autonomy Contract, Latest Stage.
- Read `FillMCTrigEffHists.cxx` (871 lines), `FitMCSinglesEffcy.cxx` (245), `plot_mc_trig_eff.cxx`
  (1853, READ ONLY — owned by another agent), `RDFBasedHistFillingData.cxx:625-655`
  (`EvaluateSingleMuonEffcyPtFitted` — the data-side ε lookup I mirror),
  `RDFBasedHistFillingPP.cxx:300-342` / `PbPb.cxx:669-712` (data fit-file paths),
  `pipelines/run_mc_trigeff_round7.sh` (driver template).
- Verified all inputs exist: nominal MC hists+fits for pp_full / overlay, both WPs (written today
  21:14–21:22); data fit files with 44 (pp) / 264 (pbpb, 6 centralities) TF1s; data 2D `_divided`
  fallbacks present incl. `_ctr0_5`.
- ROOT 6.34.04 via `bash -lc 'source /usatlas/u/yuhanguo/setup.sh'`.

### Step 1 (done) — code changes in FillMCTrigEffHists.cxx / FitMCSinglesEffcy.cxx
`FillMCTrigEffHists(sample, do_step3, use_tight_wp, do_step4, do_sanity, corrected_mc, sf_closure)`
(two NEW trailing args, both default false ⇒ nominal behaviour untouched).
- `SampleConfig` gained `data_fit_file_tmpl / data_hist_file_tmpl / data_ctr` (+`SubstWP`).
- New `DataEffEvaluator` (ε_data from the data T&P TF1s + 2D gap fallback, MCEffEvaluator's
  clamp/cap/floor guards) and `SFEvaluator` (SF = ε_data/ε_MC, pathology cap [0.01, 20] with
  counting + <SF>/σ/min/max printout).
- Steps 1/3/4 now build the numerator weight as `w·SF/ε` with `SF ≡ 1` and `ε = ε_MC` in nominal
  mode; error terms generalised to `A = Σa²`, `B = Σa²·ε_MC^nominal`, `covP = Σ2a₁a₂`,
  `covQ = Σ2a₁a₂ε_MC,1ε_MC,2` (identical to the old expressions at SF=1).
- Corrected mode skips Step 2 and the numl1/numhlt numerators (no data L1-only reference).
- `FitMCSinglesEffcy(sample, use_tight_wp, corrected_mc, sf_closure)`; object names INSIDE the
  fit file unchanged so `MCEffEvaluator` loads ε_corr with no special casing. Nominal fit-PNG
  filename left byte-identical (pre-existing quirk: nominal Medium overwrites the Tight PNG —
  NOT fixed here, out of mandate); corrected PNGs carry `_corrected[_sfclosure]<wp>`.

### Step 2 (done) — SF ≡ 1 CLOSURE TEST of the whole corrected code path (overlay, Tight)
`sf_closure=true` forces SF ≡ 1, so the corrected chain must reproduce the nominal one exactly.
Artefacts `*_corrected_sfclosure*` (nothing nominal touched). Bin-by-bin comparison
(`scratchpad/cmp_hists.C`, ALL cells incl. under/overflow):
- Step 1 : 16 hists, worst relative difference **0.000e+00** (exact)
- Step 3 : 20 hists, worst **4.62e-16** (abs 2.5e-21, `h_mc_dr_full_vs_pt_eta_errB`)
- Step 4 : 24 hists, worst **3.89e-15** (abs 6.1e-18, `h_mc_single_dr_full_vs_pt_eta_covQ`)
⇒ the refactor is exact to machine epsilon; the nominal physics is unchanged, and any
corrected-vs-original difference measured later is physics, not a code artefact.
Guard rates seen on the overlay (Tight, Step 1): 18.54% q·η-gap 2D fallbacks (identical count
for ε_MC and ε_data — same gap regions), ε_MC floor 0.24%, ε_data floor 1.95%.

### Step 3 (done) — BUG FOUND AND FIXED: `BayesDivide` is undefined for the corrected numerator
First corrected run of the overlay produced garbage forward fits (ε_corr(4 GeV) = 0.4500 in
q·η ∈ [2.0,2.2) μ⁻ where data = 0.1017). Root cause, from the fit log:
`Error in <TROOT::TEfficiency::CheckConsistency>: passed TEfficiency objects do not have
consistent bin contents` → `Error in <TGraphAsymmErrors::Divide>: passed histograms are not
consistent` → `Warning in <Fit>: Fit data is empty`.
The corrected numerator `N = Σ_fired w·SF` **exceeds** the denominator `D = Σ_all w` wherever
SF > 1, so `TGraphAsymmErrors::BayesDivide` rejects the pair and returns an EMPTY graph;
`TGraph::Fit` on an empty graph is a **no-op that still returns status 0** and leaves the TF1 at
its INITIAL parameters — 0.9/(1+e⁰)·1 = **0.4500** at pT = 4, exactly the value observed. Three+
q·η bins were affected per (sample, WP). A silent, badly wrong turn-on.
**Fix (3 parts):**
1. `FillMCTrigEffHists` corrected Step 1 now also books `errA`/`errB` for the numerator
   (`(w·SF)²` and `(w·SF)²·ε_MC`) on the pt/eta/phi 1D and the pT×q·η 2D.
2. `FitMCSinglesEffcy::CorrectedEffGraph()` builds the corrected graph explicitly:
   `eff = N/D`, `e = sqrt(A−B)/D` — the CONDITIONAL (binomial-correct) error, since
   `Var(N) = Σa²p(1−p)` with `a = w·SF`, `p = ε_MC` estimates as `A − B` from the fired muons
   (unweighted limit: `A−B = D·eff(1−eff)` ✓). k=n boundary → "1/n rule" on `n_eff = (D/e_D)²`.
   Zero-width bins (the data binning's duplicated 8.0 GeV edge) and empty denominators skipped,
   as BayesDivide does in the nominal path. NOMINAL path still uses BayesDivide, untouched.
3. Loud failure: the fitter now also counts `Npts == 0` or `ndf <= 0` as a FAILED fit, and the
   driver script aborts unless the fit stage prints "fits, 0 failed".
The plot macro's 1D corrected Step-1 series uses the same estimator + conditional error.

### Step 5 (done) — CONCURRENCY HAZARD: the pp_full single-muon NTP was being rewritten under me
The first clean pp_full corrected run (22:18–22:30) silently used a **partially written** input:
`pythia_fullsim_full_sample/muon_pairs_..._mc_trig_single_muon_full.root` was being regenerated
at that moment by the orchestrator's round-7 `run_pythia_fullsim_single_muon_mc_trig_full_sample.sh`
(adding `n_vtx`; PID 846943, finished **22:56:29**, 1.02 GB). Symptom: the closure Step-1
denominator came out **7% low at pT ≈ 4 GeV and orders of magnitude low above 20 GeV** versus the
nominal 21:17 hists — a truncated TTree, not a physics effect. Detected only because the SF ≡ 1
closure compares against the nominal file; without that test it would have propagated silently.
The pair file (Steps 3/4) was untouched (mtime 2026-07-22). The overlay single-muon NTP finished
at 21:54 and its closure was EXACT against the 21:14 nominal ⇒ the NTP rewrite changes no
kinematics, it only adds `n_vtx`.
**Action:** waited for the NTP to complete, then re-ran the ENTIRE corrected chain
(fill + closure + plots, both samples, both WPs) on the final input.
**For the orchestrator:** the nominal round-7 pp_full Step-1 hists/fits (21:17/21:18) predate this
NTP rewrite. If they are regenerated, the corrected products must be regenerated too
(`run_mc_trigeff_corrected.sh`) — they must share one input.

### Step 6 (done) — SECOND CONCURRENCY HAZARD: the orchestrator edited the SAME file mid-task
At **23:07:02** `Analysis/RDFBasedHistFilling/FillMCTrigEffHists.cxx` **in the main working tree**
(the file this task was told it OWNS) was modified by the orchestrator, recompiled 23:08, and the
NOMINAL artefacts regenerated 23:12–23:16. The change is the round-7 **FORWARD LOW-pT VETO**
(contract item 3): `kVetoFwdLowPt = true`, a muon is rejected iff `pT < 7 GeV AND q·η < −2`,
applied to **Steps 2/3/4 only** (`kFwdVetoLeg` in `sel_pair_legs`, `kFwdVetoPair` in the Step-3
selection); Step 1 and the sanity mode stay unvetoed.
Detected by the SF ≡ 1 closure: the corrected Step-3 denominator came out ~9% ABOVE the nominal
one (`h_mc_dr_zoom_vs_pair_pt_denom` integral 45.97 vs 42.13) — my copy still lacked the veto.
Without the closure test this would have compared a vetoed nominal against an unvetoed corrected
sample and reported the difference as physics.
**Action:** re-synced the veto into the worktree copy VERBATIM (constants + both selection
strings + the header block, plus master's `weighted%%`→`weighted%` printf fix), verified that the
only remaining differences vs master are the intended corrected-MC additions, and re-ran the whole
chain. **The corrected products must always be regenerated whenever the nominal selection changes
— that is exactly what `run_mc_trigeff_corrected.sh` is for.**

### Deliverables written (paths)
Code (all in the worktree; cherry-pick, do NOT merge the branch):
- `Analysis/RDFBasedHistFilling/FillMCTrigEffHists.cxx` (edited)
- `Analysis/RDFBasedHistFilling/FitMCSinglesEffcy.cxx` (edited)
- `Analysis/plotting_codes/trig_effcy/mc_based/plot_mc_trig_eff_corrected.cxx` (NEW)
- `Analysis/pipelines/run_mc_trigeff_corrected.sh` (NEW)
ROOT artefacts, per sample dir (`pythia_fullsim_full_sample/`,
`pythia_fullsim_hijing_overlay_test_sample/`), per WP (`` / `_medium_wp`):
- `mc_trig_eff_hists_<label><wp>_corrected[_step3|_step4].root`
- `single_mu_effcy_pT_fit_mc_corrected<wp>.root`
- closure twins with `_corrected_sfclosure`
Plots, under `<plot base>/step{1,3,4}_corrected_mc/[medium/]`, plot base =
`dimuon_data/plots/pp_trigger_efficiency/mc_based/` (pp_full) and
`.../pbpb_trigger_efficiency/mc_based/` (overlay).

### Step 4 (done) — SECOND BUG in the new estimator: eff = 0 +- 0 destroys the fit
First run with the new `CorrectedEffGraph` gave `pp_full` Tight, q·η ∈ [1.3,1.6) μ⁺:
**chi2/ndf = 30626.9/30, eps(30 GeV) = 1.396** (unphysical). Cause: a high-pT bin with
denominator weight > 0 but NO fired muons (A = B = 0) fell into the boundary branch and got
`e = (N/D)/n_eff = 0` → an `eff = 0 ± 0` point is infinitely constraining. Fix: the "1/n rule"
must cover BOTH binomial boundaries, `e = max(eff, 1)/n_eff` (k = n → e ≈ 1/n; k = 0 → e ≈ 1/n).
Applied in `FitMCSinglesEffcy::CorrectedEffGraph` and in the plot macro's twin.
(The nominal path never saw this: `BayesDivide` gives a k = 0 point a proper Bayesian upper error.)
Also fixed in the driver: `root ... | tee /dev/stderr` re-opens the log through
`/proc/self/fd/2` at offset 0 and OVERWRITES it — the fit output is now written to its own
`fit_<sample>_<wp>*.log` and then echoed + grepped.
Also fixed in the plot macro: `Form()` was called with a `%.0s` conversion against a `double`
argument in the pull-pad off-scale label (undefined behaviour); the format now has two
conversions and the extra argument is simply unused.
