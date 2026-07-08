# Upfront hadronic/fake background reduction (working points, |d0|, Δp/p)

**Mode:** Implementation. **Created:** 2026-07-06. **Session:** "template fit — upfront bkg reduction".
**Parent:** `low_mass_dimuon_template_fit.md` (sub-doc **E** of the factorized template-fit set).
**Siblings:** `tf_bkg_composition_normalization.md` (C — provenance classifier, tight-vs-medium
study, data/MC), `tf_dpop_fake_muon_program.md` (D — the Δp/p *fit-and-subtract* yield program;
THIS doc is the *cut-upfront* counterpart), `analysis_overview.md` §6 (three low-mass backgrounds).
**Reviewer rules:** ParamsSet/macro C++ → `/review-analysis-code`; plots → `/review-plot`.

---

## Objective
Assess whether the enhanced **hadronic (π/K decay-in-flight, punch-through) and fake** muon
background in the **low pair-mass region** can be reduced **upfront by tighter muon selection**,
rather than relying on the pair-level template fit to extract and subtract it. Three
complementary handles, pursued in parallel:
1. **Working points (WP).** Document exactly what the **Medium** and **Tight** muon quality
   working points require (from the ATLAS/Athena source), so we can judge whether an
   *additional* cut (e.g. a stricter track-fit χ²/normalized-χ²) on top of the WP would remove
   hadronic/fake muons.
2. **|d0| discrimination.** Real muons (open-HF dominated) have a **displaced** production point
   (b/c lifetime) → finite |d0|; hadronic/fake muons are **prompt on the ID scale** → pile up at
   very small |d0|. Plot |d0| by provenance (log-x, small-|d0| focus) to see whether a **lower**
   |d0| cut is a useful discriminant.
3. **Δp/p distribution.** Characterize Δp/p by provenance (all single-muon cuts *except* Δp/p),
   to see how well the nominal `Δp/p<0.12` cut separates hadronic/fake from real, and its **muon-pT
   dependence** (the boost argument, below).

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation
The low-mass OS spectrum sits on three backgrounds (`analysis_overview.md` §6): uncorrelated
combinatoric (charge-symmetric), correlated gluon-splitting/open-HF G (real μμ), and **fake/hadronic**
μ (≥1 π/K decay-in-flight, punch-through, or mis-associated fake). The template-fit program
(sibling docs A/D) removes these downstream. This doc asks the orthogonal question: **how much can
be removed at the muon-selection stage**, which is cleaner (no fit model, no double-subtraction
risk — see sibling D "double-subtraction") and directly protects the pT spectrum → R_AA(pT) from a
pT-dependent contamination.

### 2. |d0| physics (why it discriminates, and the pT dependence)
`muon_d0` = transverse impact parameter (closest transverse approach of the extrapolated muon
track to the primary vertex/beamline), **unit mm**. Nominal analysis applies an *upper* bound
`|d0| < d0cut = 2 mm` (`ParamsSet.h:365`); this study evaluates a possible *lower* bound.
- **Real (open-HF).** The b/c hadron is produced at the PV but travels before decaying —
  cτ ≈ 460 µm (B), ~120–300 µm (D), × Lorentz boost βγ → mm-scale transverse flight L_xy. The muon
  is born at that **displaced secondary vertex** → track does not point back to the PV → a
  characteristically **finite |d0|**. "Bounded below" = a statistical scale set by the decay
  length (NOT a hard floor: a muon emitted along the flight direction, or from an early decay, can
  have small d0), so the *ensemble* has its weight pushed to finite |d0| and does not pile up at 0.
- **Hadronic (π/K decay-in-flight, punch-through).** π±/K± originate essentially at the PV (prompt
  on the ID scale); the decay-in-flight muon track largely still points at the beamline (small
  kink for energetic parents), punch-through hadrons are prompt → |d0| set by **tracking
  resolution, not a decay length** → **peaks at very small |d0|**.
- **Fake (no truth match).** Random ID-track↔MS-segment associations of prompt tracks → likewise
  resolution-dominated → small |d0|.
⇒ Small-|d0| region is dominated by hadronic+fake; real-HF is suppressed there. A **lower |d0|
cut** would preferentially remove hadronic/fake while keeping most real-HF — orthogonal to Δp/p.

**pT dependence (ties |d0| to Δp/p — 3 effects to display):**
(a) **Hadronic survival of the Δp/p cut rises with pT.** Δp/p=(p_ID−p_MS)/p_ID flags the ID-vs-MS
    imbalance from a decay kink; π/K decay-in-flight → large positive Δp/p → removed by
    `Δp/p<0.12`. But in-flight decay probability ∝ 1/(γcτ): higher muon pT → larger parent γ →
    longer lab lifetime γτ → decays **later/deeper** (past the ID, in/beyond the calorimeter) →
    smaller kink → **smaller Δp/p → more likely to survive**. So surviving hadronic contamination
    grows with pT. The pT-binned Δp/p THStack is the direct diagnostic.
(b) **|d0| is complementary exactly there.** |d0| keys on the presence/absence of a decay-length
    displacement, independent of the kink → can still remove high-pT hadronic muons that evade Δp/p,
    provided real-HF keeps its displacement at that pT.
(c) **Does d0 itself depend on pT?** Real-HF d0 is set by transverse decay kinematics more than by
    boost → only *weakly* pT-dependent (why impact-parameter HF taggers work across pT). Hadronic/fake
    d0 is resolution-dominated, and tracking d0 resolution *improves* with pT (straighter tracks) →
    their small-|d0| peak *narrows* at high pT. If real-vs-hadronic |d0| separation shifts with pT →
    study finer pT bins and/or a **pT-dependent |d0| cut**.
**Guardrail:** an unaddressed pT-dependent hadronic contamination would distort the measured pT
spectrum → bias R_AA(pT). The pT-binned THStack is what catches it.

### 3. Provenance classification (4-way; from truth, reco-seeded)
Per reco muon (branches verified present in the NTUP):
- **fake** = `muon_truth_prob ≤ 0.5` (no truth match).
- **hadronic** = `prob>0.5 & (|muon_truth_id|≠13` [punch-through hadron] `OR (|id|==13 &
  muon_truth_IsPrimary==0)` [decay-in-flight → Geant secondary muon]`)`.
- **real** = `prob>0.5 & |muon_truth_id|==13 & IsPrimary==1`. **Split real into (per user):**
  - **real-HF** — parent is a c- or b-hadron (displaced).
  - **real-prompt** — parent is W/Z/γ*/H, τ, or quarkonium (J/ψ, ψ′, Υ, …) (d0≈0).
  HF-vs-prompt via truth navigation (NO MCTruthClassifier origin branch exists): reco muon i →
  `j=muon_truth_index[i]` (skip if <0) → parents `truth_parents[j]` (parent **barcodes**) → map
  each barcode to its index via `truth_barcode` → parent pdgId `truth_id[k]` → classify:
  HF if the heaviest quark of |parentId| is c(4) or b(5) (mesons: `(|id|/100)%10`; baryons:
  `(|id|/1000)%10` ∈ {4,5}); prompt if |parentId| ∈ {23,24,22,25,15} or a quarkonium
  (443,100443,553,100553,…). **VALIDATION REQUIRED:** confirm real-HF shows a displaced truth
  production vertex (transverse `sqrt(truth_vtx_x²+truth_vtx_y²)` ≫ 0 or non-zero d0) while
  real-prompt sits at ≈0 — this cross-checks the parent-pdgId logic. Empirically determine whether
  `truth_parents` holds barcodes or indices (barcodes expected) before trusting the map.
- Labels on plots: **"real HF"**, **"real prompt"**, **"hadronic"**, **"fake"** — NEVER "prompt"
  for the HF muons (memory `feedback_hf_muons_not_prompt`).

### 4. Muon selection (reco-seeded, applied in-macro on the raw NTUP)
Medium WP = `muon_quality & (1|8|32|256) == (1|8|32|256)` (combined|Medium|IDCuts|MuonCuts);
Tight adds bit 16 → `(1|16|32|256)`. Kinematic: `pt>4 GeV`, `|η|<2.4`. Analysis cuts:
`Δp/p<0.12` (`deltaP_overP_thrsh`), `|d0|<2 mm` (`d0cut`), `|z0 sinθ|<2 mm` (`z0cut`).
- **|d0| plot (task 2):** apply **ALL** muon cuts incl. Δp/p and the |d0|<2mm upper bound, then
  histogram |d0| (log-x, small-|d0| focus, e.g. [1e-3, 2] mm).
- **Δp/p plot (task 3):** apply all single-muon cuts **EXCEPT Δp/p** (keep d0/z0/quality/pt/η),
  then histogram Δp/p. **No ntuple reprocessing needed** — the raw NTUP retains all muons and the
  Δp/p cut is analysis-level (existing `dpop_templates.C` already plots Δp/p with no Δp/p cut).

### 5. Plot specification (both |d0| and Δp/p)
Two plotting **modes**, each in **pT-integrated** and **coarse single-muon-pT-binned** versions
(4 plots per variable per sample; × {|d0|, Δp/p} × {signal fullsim, HIJING overlay}):
- **Mode 1 — unity-normalized overlay:** each provenance class normalized to unit area, 4 lines
  overlaid (real-HF, real-prompt, hadronic, fake). Shows the *shape* discrimination.
- **Mode 2 — THStack, differential cross-section:** stacked dσ by provenance (pT-hat-slice
  AMI-weighted, nb→pb per `fill_weighted_fullsim.C`; use `Scale(NORM,"width")` per memory
  `feedback_differential_scale`). Shows the *absolute* composition sitting at small |d0| / surviving
  Δp/p, and its pT dependence.
- **Axes:** |d0| log-x (distribution → also log-y per memory `feedback_log_scale_plots`); Δp/p
  linear-x over ~[−0.6,0.6] (log-y for the distribution). Subplot grid per memory
  `feedback_subplot_layout` (N≤3 → 1 row).

### 6. Single-muon pT binning (new, in ParamsSet.h)
Add a **coarse single-muon pT binning** to `ParamsSet.h` as the single source of truth (alongside
`pair_pt_coarse_bins`), guideline ≈ half the pair-pT scale (pairs 8/15/27/50/150 → single-muon
~4/7.5/13.5/25/75), respecting `pt>4`. **The subagent optimizes the edges** from the actual reco
muon-pT spectrum (signal+overlay) to ensure coverage of the probable pT range, sensitivity to a
real pT dependence, and adequate per-bin stats; then registers them in ParamsSet.h with a comment
that this is the source of truth (mirror `project_coarse_pair_pt_binning`). Both the |d0| and Δp/p
macros read these edges (source-of-truth; do not re-invent per plot).

### 7. Negative constraints
- Do NOT reprocess ntuples for the Δp/p "no-cut" sample — it already exists in the raw NTUP.
- Do NOT label real-HF muons "prompt".
- THStack must use the SAME AMI/pT-hat weighting as the validated `fill_weighted_fullsim.C`
  (nb→pb); do not re-derive a different normalization.
- |d0| and Δp/p macros/plots live in SEPARATE directories under the template-fit plots area; they
  are parallelizable EXCEPT for the shared ParamsSet single-muon pT binning (do that first).
- The single-muon pT binning is binned in **single-muon** pT (≈ half pair-pT), NOT pair pT.

## Context
- Samples: `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/Pythia_5p36TeV_{pp,pn,np,nn}_hQCD_DiMu_pTH*.FullSimPP24.NTUP.root`
  (24 slices, pp-condition signal) and `.../pythia_fullsim_hijing_overlay_test_sample/*.FullSimHIJINGOverlayPP24.NTUP.root`
  (6 slices, overlay). Tree `HeavyIonD3PD`.
- Branches (verified): `muon_pt` (MeV, signed=charge), `muon_eta`, `muon_phi`, `muon_quality`,
  `muon_deltaP_overP`, `muon_d0` (mm), `muon_z0`, `muon_truth_prob`, `muon_truth_id`,
  `muon_truth_IsPrimary`, `muon_truth_index`; truth collection `truth_id`, `truth_barcode`,
  `truth_parents` (vector<vector<int>>), `truth_vtx_{x,y,z}`.
- Reusable machinery: `dimuon_data/plots/template_fitting/dpop_fake_muon_20260624/code/dpop_templates.C`
  (reco-seeded Δp/p templates, 3-way classify — extend to 4-way + pT bins + d0/z0 cut);
  `bkg_mc_provenance_20260624/code/fill_weighted_fullsim.C` (AMI pT-hat weighting nb→pb).
- Athena WP source: `cd ../SkimCode && source setup_25.sh && which athena`; MuonSelectorTools
  quality definitions in `/cvmfs/atlas.cern.ch/repo/sw/software/25.2/AthAnalysis/25.2.89/InstallArea/x86_64-el9-gcc14-opt/src/`.
- Plot output area: `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/template_fitting/`.

## Scope
In: (1) `muon_working_points.md` repo doc + χ² recommendation; (2) single-muon pT binning in
ParamsSet.h; (3) |d0| plots (2 modes × pT-int/binned × 2 samples) in a `d0_discrimination_<date>/`
dir; (4) Δp/p plots (same matrix) in a `dpop_distribution_<date>/` dir. Out: actually *changing*
the analysis selection (this is a study to decide whether to); the Δp/p yield fit (sibling D).

## Design Decisions
- **4-way provenance split (real-HF / real-prompt / hadronic / fake)** — user-approved 2026-07-06
  (over the 3-way), so the "displaced" feature is cleanly attributable to HF and a prompt/quarkonium
  small-|d0| sub-population cannot masquerade as discrimination-limiting. Adds truth-parent
  navigation (§3) with a production-vertex validation.
- **No ntuple reprocessing** — the raw NTUP retains pre-Δp/p-cut muons (verified); the standalone
  reco-seeded macro pattern (read NTUP, cut in-macro) is reused. Resolves the user's "may need a
  no-Δp/p single-muon-tree mode" concern.
- **Single-muon pT binning determined by the subagent from the spectrum**, saved to ParamsSet.h as
  source of truth (per user; mirrors `project_coarse_pair_pt_binning`).

## Implementation Plan
1. **Doc + Physics Procedure** (this file). DONE.
2. **T1 — WP repo doc (task 1)** — Athena source dig → `Analysis/docs/references/muon_working_points.md`
   (Medium vs Tight requirements, esp. any χ²/normalized-χ² and hit requirements) + a recommendation
   on whether a stricter χ² on top of the WP would help. → `/review-analysis-code` (doc/physics check,
   quote §1). INDEPENDENT, parallel. **DONE** (commit 80de5a4).
3. **T2 — single-muon pT binning in ParamsSet.h (task 2/3 prereq, §6)** — determine edges from the
   muon-pT spectrum, add to ParamsSet.h, ACLiC-verify. → `/review-analysis-code` (quote §6).
   Must precede T3/T4 (shared file). **DONE** (commit cec2416).
4. **T3 — |d0| plots (task 2, §2/§4/§5)** — new dir `d0_discrimination_<date>/`; 4-way classify,
   both modes, pT-int + binned, signal + overlay. → `/review-plot` (quote §2,§3,§5). **DONE** (PASS iter 2).
5. **T4 — Δp/p plots (task 3, §4/§5)** — new dir `dpop_distribution_<date>/`; same matrix. →
   `/review-plot` (quote §2c,§3,§5). Parallel with T3 after T2. **DONE** (PASS iter 2).

## Progress Log
- 2026-07-06 — Step 1: doc created (factorized sub-doc E). Physics Procedure established: |d0| &
  Δp/p discrimination + pT dependence (§2), 4-way provenance with truth-parent HF/prompt split
  (§3), reco-seeded no-reprocessing approach (§4), 2-mode × pT-int/binned × 2-sample plot matrix
  (§5), single-muon pT binning to ParamsSet.h (§6). Verified in code: `muon_d0`/`muon_deltaP_overP`
  branches, cut values (Δp/p<0.12, |d0|<2mm, medium=1|8|32|256/tight+16), truth-nav branch types,
  samples (24 signal + 6 overlay slices), reusable macros.

- 2026-07-06 — **T2 (single-muon pT binning) DONE** (subagent `_sub_binning_1`, merged+deleted).
  Measured the reco muon-pT spectrum (medium sel) per provenance in signal + overlay → coarse
  edges **{4,8,14,25,100} GeV** (4 bins, ~½ the pair-pT scale, open-ended top). Added to
  `ParamsSet.h` (`single_mu_pt_coarse_bins` + `N_COARSE_SINGLE_MU_PT_BINS=4`), matching the
  `pair_pt_coarse_bins` idiom; ACLiC-verified (test prints correct edges). Committed `cec2416`.
  Raw per-bin counts (overlay all/real/had/fake): 4–8 57167/52018/4723/426; 8–14 19324/18602/682/40;
  14–25 8850/8582/236/32; 25–100 4506/4438/63/5. Hadronic healthy every bin; fakes thin only in the
  top bin (expected steep fall; test-subset raw counts, scale up with full production + weighting).

- 2026-07-06 — **T1 (working points, task 1) DONE** (subagent `_sub_wp_doc_1`, merged+deleted).
  Created repo doc `Analysis/docs/references/muon_working_points.md` from AthAnalysis 25.2.89
  `MuonSelectionTool.cxx` (cited file:lines). Findings: Medium & Tight both require combined type,
  IDCuts hit reqs, ≥2 precision layers, flat |q/p significance|<7. **Tight adds (Medium lacks):**
  (i) normalized track-fit **χ²<8** (`reducedChi2`), (ii) pT/η-binned tight cuts on ρ′ (ID/MS
  momentum balance) + q/p significance, (iii) drops Medium's single-layer central exception. The
  bitmask is packed in the SKIM (`SkimCode/.../TrigRates.cxx:1193-1203`), not Analysis.
  **χ² recommendation:** Medium does ~no χ² selection → a stricter `reducedChi2<8` and/or
  tightening our existing flat `|Δp/p|<0.12` (same physics as ρ′) is non-redundant and should cut
  additional hadronic/fake; but the dominant low-mass bkg is real HF muons (unaffected), so gains
  target only the fake/decay-in-flight sub-component, and Tight-style cuts cost pT-dependent
  efficiency (worst at low pT = our signal region) → would force re-deriving the reco-eff
  placeholder. Path: add cuts one variable at a time, validate purity vs real-μ efficiency loss
  against Pythia-fullsim truth. NOT git-committed yet (bundle with the factorization docs commit).

- 2026-07-06 — **T3 (|d0|, task 2) + T4 (Δp/p, task 3) plots BUILT & VALIDATED** (subagents
  `_sub_d0_plots_1`, `_sub_dpop_plots_1`; reco-seeded standalone macros; 8 PNGs each; SUMMARY.md in
  each dir). Dirs `dimuon_data/plots/template_fitting/{d0_discrimination_20260706,dpop_distribution_20260706}/`.
  Both: 4-way provenance (real-HF/real-prompt/hadronic/fake) via truth-parent BFS chain nav, both
  modes (unity-norm overlay + AMI-weighted nb→pb THStack), pT-int + coarse single-μ-pT-binned,
  signal + overlay. **/review-plot = the remaining gate** (see Latest Stage).
  - **CLASSIFICATION BUG (|d0|) found + fixed by orchestrator:** the fill capped the truth-particle
    count by `min(sizes)` INCLUDING `truth_vtx_*` (a separate ~30%-shorter VERTEX collection) →
    starved the parent barcode→index map → real-HF muons collapsed into the real-prompt residual.
    Fix: count over `truth_id/barcode/parents` only; upgraded BOTH macros to a full BFS ancestry
    chain-walk (any c/b-hadron ancestor ⇒ real-HF; catches b→c→μ and B→J/ψ→μ). Also optimized the
    Δp/p macro's per-event barcode map to lazy/unordered (identical output, ~5× faster). GOTCHAS
    (recorded in the SUMMARYs): `truth_parents` needs a runtime `GenerateDictionary` or reads size 0;
    `truth_vtx_*` is not particle-indexed → displacement validated via reco |d0|.
- 2026-07-06 — **Doc E FACTORIZATION committed** (`522484e`) + WP doc (`80de5a4`) + single-μ binning
  (`cec2416`). Plot dirs are data-area (not git-tracked).

- 2026-07-07 — **CERTIFIED presentation-ready** (`/review-plot` PASS, log
  `review-plot-20260707-120911-upfront-d0-dpop-certify.md`; 0 CRITICAL/0 WARNING; C1–C6 incl. C5
  provenance + C6 MC-weight). Round of user-driven corrections, all applied to BOTH macros +
  regenerated + re-reviewed:
  1. **Overlay classification — INDEX gate.** The generator-block *barcode* cutoff (for the parent
     map) recovered overlay real-HF 5%→~90%, but ~1724 overlay muons with `muon_truth_index≥lim`
     (matched to the HIJING/Geant4 block; stored `muon_truth_barcode≤200000` but index out of block —
     an overlay index/barcode inconsistency) were still mis-labelled "real-prompt residual". Fix: a
     real muon now REQUIRES `muon_truth_index∈[0,lim)` (matched INTO the Pythia-signal block); else →
     hadronic/background. Overlay real now 91.7% HF, **0 unresolved**, real-prompt 6167 (med |d0| 0.017,
     smallest), hadronic 2390→4114. Gate is in the COMMON path (applied to pp24 too; a no-op there —
     pp24 real counts unchanged, 0 unresolved). [[reference_muon_truth_provenance_nav]] updated.
  2. **Per-event MC weights on BOTH modes.** dpop mode-1 (unit-normalized) had used RAW counts →
     fixed to the weighted histos (d0 mode-1 already weighted). Violated the standing rule (memory
     `feedback_hf_muons_not_prompt`); enforcement added as **criterion C6** in
     `conventions/physics-results-review.md` (read by /review-plot AND /review-analysis-code),
     committed `b2e983d`.
  3. **Recolor** red/blue/kGreen+2/magenta; **combined 3-way** variant (`combined/` subdir);
     **descriptive sample labels/filenames** `fullsim_pp24` / `hijing_overlay` (on-plot "pythia
     fullsim, pp24" / "pythia fullsim, hijing overlay, Pb+Pb23"); **Δp/p cut-line label** fixed
     (reads "Δp/p = 0.12" at the top strip, not clipped/buried).
  Non-blocking INFO left as-is (legible): d0 top subtitle grazes the top frame; dpop 2×2 sample-label
  leading char grazes the pad-left margin. **FOLLOW-UP:** `bkg_mc_provenance_20260624` uses the SAME
  classifier and lacks the index gate → its overlay real/prompt is inflated (HIJING UE counted as
  real); apply the same index gate there if its overlay composition is used quantitatively.

- 2026-07-07 — **[SUPERSEDED — the index gate was a MISTAKE, REVERTED same day; see the INDEX-GATE
  REVERT entry below. This entry is kept as history only.]** **FOLLOW-UP DONE: index gate applied to
  `bkg_mc_provenance_20260624`** (subagent
  `_sub_bkgprov_indexgate_1`, merged+deleted). The same Pythia-signal-block **index gate** was added
  to BOTH classifiers there — `bkg_mc_provenance.C` and `fill_weighted_fullsim.C` (read
  `muon_truth_index`/`muon_truth_barcode`/`truth_barcode`; per-event `lim`=first `truth_barcode>200000`;
  gate `ti<0 || ti>=lim || tbc>200000 → hadronic`, placed AFTER the fake test, BEFORE the prompt test).
  Backup first → `backup_pre_indexgate_20260707/` (plots 6, plots_tight 12, 8 .root, code). Recompiled,
  reran all fills + regenerated all plots (medium `plots/` now 12 PNG — full 3-plotter regen added the
  make_plots/weighted-meeting outputs that pre-fix had only 6 of; `plots_tight/` 12 PNG).
  **No-op / shift verification (backup vs new histogram integrals):**
  - **SIGNAL (pp24), all 4 files (prov med/tight, wfill med/tight): EVERY histogram Δ = 0 — EXACT no-op.**
  - OVERLAY: real drops / hadronic rises by identical amount; fake & totals UNCHANGED (gate is after the
    fake test): prov MED op prompt −50/had +50, ss −49/+49; prov TGT op −45/+45, ss −47/+47; wfill MED op
    real −11.15/had +11.15, ss −13.67/+13.67; wfill TGT op −10.40/+10.40, ss −13.26/+13.26. Single-muon
    overlay medium: prompt 96.57%→94.37%, hadronic 3.05%→5.25% (1724 real→had); pairs-overlay totals
    unchanged (19903). A separate read-only agent is diffing backup-vs-new PLOTS (esp. `data_vs_fullsim_*`,
    which must stay identical — no HIJING there) to confirm no adverse effect on the data-vs-fullsim view.
  Sibling `tf_bkg_composition_normalization.md` (C) updated to note its overlay real/prompt numbers were
  inflated by the same missing gate (now fixed). [[reference_muon_truth_provenance_nav]].

- 2026-07-07 — **INDEX-GATE REVERT, step 1 (investigation) DONE.** Diagnostic
  (`scratchpad/diag_untraced.C`, overlay, mirrors d0 selection+lim EXACTLY, correct real def no gate):
  fake=299 hadronic=2390 **real=75659** (traced idx∈[0,lim)=73935; **untraced idx<0||≥lim=1724 =
  2.28% of real**; idx<0=3, idx≥lim=1721). ALL 1724 untraced reals: `truth_id=±13`, `IsPrimary=1`,
  `muon_truth_barcode≤200000` (1473 ≤47k, 251 (47k,200k], **0 >200k** → none Geant4). Their
  `muon_truth_index` is 8k–69k ≫ `lim`(~400–1400): the truth array is NOT barcode-sorted; `[0,lim)`
  captures the embedded **Pythia-signal** block, and HIJING underlying-event **generator primary
  muons** (bc 1–47k) sit at high array positions → untraceable within `[0,lim)`. 58/1724 have
  `bcInGenBlock=1` = barcode-duplication cases (why the `[0,lim)` parent-map cutoff is mandatory).
  **⇒ CONFIRMED: the untraced reals are HIJING primary muons (real), NOT Geant4/hadronic** — validates
  the user's correction and the label **"real (untraced HIJING)"**. Signal (pp24) untraced≈0 (no HIJING
  UE dumping high-index generator muons; consistent with the earlier EXACT signal no-op).

- 2026-07-07 — **INDEX-GATE REVERT IMPLEMENTED + regenerated + /review-analysis-code PASS**
  (log `review-analysis-code-20260707-202032-index-gate-revert-5class.md`, iter 1, 0 CRIT/0 WARN;
  C5 provenance + C6 weights verified). Changes:
  - **d0_discrimination.C + dpop_dist.C(+plot_dpop_dist.C):** removed the index-gate line; restored the
    correct hadronic def (`prob>0.5 & (|id|≠13 || IsPrimary==0)`); KEPT the `[0,lim)` parent-map cutoff
    (HF/prompt sub-split only); added a 5th overlay-only class **"real (untraced HIJING)"** (real muons
    with `muon_truth_index<0||≥lim`), color kOrange+1. Reordered classes 0=realHF 1=realprompt
    2=realuntracedHIJING 3=hadronic 4=fake (reals contiguous at stack bottom); combined 3-way "real (all)"
    = 0+1+2.
  - **bkg_mc_provenance.C + fill_weighted_fullsim.C:** restored VERBATIM from `backup_pre_indexgate_20260707/`
    (3-way real/fake/hadronic, no index gate).
  - **Regenerated outputs (all confirmed):** d0 16 PNG, dpop 16 PNG, bkg 12+12 PNG.
    - **d0 overlay:** real HF 67768 (med|d0| 0.088) / real prompt 6167 (0.017) / **real (untraced HIJING)
      1724 (0.057)** / hadronic **2390** (0.021, smallest) / fake 299. real total 75659 = 67768+6167+1724.
      Hadronic back to 2390 (the index gate had wrongly inflated it to 4114). **Signal (pp24): 0 untraced
      HIJING** (no-op).
    - **dpop overlay** (no Δp/p cut → looser): real (untraced HIJING) pT-int N=2031 (>d0's 1724 as expected),
      med|d0| 0.084, Δp/p mean +0.035 (between real +0.014 and hadronic +0.09); signal 0.
    - **bkg overlay back to ORIGINAL:** prompt(=real) **75659 (96.57%)** / fake 299 (0.38%) / hadronic
      **2390 (3.05%)** — HIJING muons counted as real again (was 94.37%/5.25% under the gate). Signal
      prompt 362217 (99.32%) / hadronic 2408 (0.66%) unchanged. bkg "prompt" total 75659 == d0 real total ✓.
  - **Docs/memory corrected:** umbrella doc §4 (HIJING-muons-are-real general rule), sub-doc C note revised
    (index gate was a mistake, reverted; overlay numbers stand), memory `reference_muon_truth_provenance_nav`
    rewritten. Physics: the untraced-HIJING |d0| (0.057) sits between prompt (PV) and HF (displaced) → a mix
    of HIJING open-HF + prompt, as expected for real UE muons.

## Results & Observations

### |d0| discrimination (T3) — real-HF displaced, hadronic/fake at small |d0| (hypothesis confirmed)
Validation (unweighted, median |d0| [mm]; overlay AFTER the 2026-07-07 barcode-cutoff fix):
**signal** real-HF 0.104 (displaced) / real-prompt 0.030 (at PV) / hadronic 0.045 / fake 0.081 —
92% real-HF; **overlay** hadronic **0.021 (smallest)** / real-prompt 0.020 (at PV) / fake 0.039 /
real-HF 0.088 (displaced) — 90% real-HF (7.6% residual, 2.3% unresolved). Overlay composition now
MATCHES signal (HIJING adds hadronic/fake but does NOT kill HF or make prompt). Readout: (1) signal
& overlay — real-prompt peaks sharply at small |d0|, hadronic smallest, real-HF displaced (§2
confirmed); (2) a lower-|d0| cut removes hadronic/fake/prompt preferentially, orthogonal to Δp/p;
(3) real-prompt small-|d0| peak sharpens with pT (resolution) while real-HF stays broad → separation
holds at high pT. **ROOT-CAUSE (corrected):** the earlier overlay "real-HF 5%, 92% real-prompt
residual" was a barcode-duplication artifact (map built over the full merged truth record without
the generator-block cutoff), NOT physics — fixed by restricting map+trace to [0,lim), lim=first
truth_barcode>200000 (overlay HF 5%→90%). See Latest Stage + [[reference_muon_truth_provenance_nav]].

### Δp/p distribution (T4) — Δp/p is a PARTIAL discriminant; motivates |d0| as complement
Validation (signal, pT-int): real-HF mean Δp/p +0.014 (symmetric, |d0| 0.206 displaced) / real-prompt
+0.014 (|d0| 0.036 at PV) / hadronic **+0.092 (shifted positive)** / fake +0.021. Hadronic survival
of the one-sided Δp/p<0.12 cut vs single-μ pT: signal 0.565/0.569/0.603/0.530, overlay
0.506/0.579/0.614/0.516 across {4-8,8-14,14-25,25-100} GeV. Readout: real muons symmetric at ~0,
hadronic broad + shifted positive (§2c). ~50–60% of hadronic SURVIVES the cut at ALL pT (only a mild
boost-consistent rise 4-8→14-25 GeV, non-monotonic, top bin low-stat) while ~94% of real-HF is kept
⇒ Δp/p is a real but PARTIAL discriminant; a large hadronic contamination evades it at every pT →
motivates |d0| as an orthogonal complement (T3).

### Combined verdict
|d0| and Δp/p are orthogonal handles on the hadronic/fake background. Δp/p (already applied) keeps
~94% signal but leaves ~50-60% hadronic; |d0| catches the surviving hadronic (smallest |d0| in the
realistic overlay). Working points (T1): Medium does ~no χ² selection; adding `reducedChi2<8` / a
tighter Δp/p, one variable at a time, is the recommended upfront-tightening path. All three are
**complementary secondary handles**, not strong standalone cuts — the low-mass hadronic/fake yield is
modest and the dominant low-mass background is real HF (untouched by these cuts). Any adoption of a
lower/pT-dependent |d0| cut or stricter χ² is a signal-selection change → `signal_selection_change_impact.md`.

## Remaining Work
- **User decision (the deliverables now inform this):** adopt a lower / pT-dependent |d0| cut,
  a stricter track-fit χ² (reducedChi2<8), and/or a tighter Δp/p? Evidence: all three are
  complementary secondary handles (real-HF displaced vs hadronic-at-small-|d0|; Δp/p leaves
  ~50–60% hadronic; Medium does ~no χ²), NOT strong standalone cuts; dominant low-mass bkg is real
  HF (untouched). Any adoption is a signal-selection change → `signal_selection_change_impact.md`.

## Latest Stage
**✅ 2026-07-07 — INDEX-GATE REVERT COMPLETE & CERTIFIED.** Both reviewer loops PASS (0 CRIT/0 WARN):
`/review-analysis-code` (log `review-analysis-code-20260707-202032-index-gate-revert-5class.md`) +
`/review-plot` (log `review-plot-20260707-203156-index-gate-revert-5class.md`). The wrong index gate is
removed from all classifiers; correct definition restored (HIJING muons ARE real); new overlay-only class
**"real (untraced HIJING)"** added (orange) in d0/dpop; bkg classifiers restored from backup. Regenerated:
d0 16 PNG, dpop 16 PNG, bkg 12+12 PNG. Overlay: real(untraced HIJING)=1724 (d0) / 2031 (dpop, no Δp/p cut),
hadronic back to 2390 (3.05%), real total 75659 (==bkg prompt). Signal pp24: 0 untraced (no-op). Docs +
memory `reference_muon_truth_provenance_nav` + umbrella §4 (HIJING-real general rule) corrected. Plot dirs
are data-area (not git). **Remaining:** only the user's selection-tightening decision (see Remaining Work) —
the study deliverables are final. **Clear Latest Stage after commit.**

Reference (the reverted mistake, for history):
**🔴 2026-07-07 — REOPENED: the INDEX GATE was a PHYSICS MISTAKE; reverting.** User correction:
**HIJING muons ARE real muons** (HIJING simulates hard scatterings too → single-B signal, open-HF,
prompt light-meson/quarkonium/QED, and combinatoric — like Pythia). The **original** definition
`hadronic = prob>0.5 & (|id|≠13 || IsPrimary==0)`, `real = prob>0.5 & |id|==13 & IsPrimary==1` is
CORRECT; requiring "real" to be matched into the Pythia-signal block (the index gate) is WRONG.
Reason we use Pythia-only truth for reco-eff/detector-response and for the open-HF template
backgrounds is **event weights** (Pythia-AMI weights don't apply to HIJING; mixing generators for HF
would need a full truth-history rewrite), NOT because HIJING muons aren't real. This is a general
rule for HIJING/PbPb-condition work (documented in the umbrella doc).

**Two distinct barcode mechanisms — separate them:**
- **① Index gate (`m_ti<0 || m_ti>=lim || m_tbc>200000` → hadronic): the MISTAKE → REMOVE.**
  `d0_discrimination.C:186`, `dpop_dist.C:181`, and inside `classify()` in `bkg_mc_provenance.C` +
  `fill_weighted_fullsim.C`. Removing it restores the correct real/hadronic/fake axis (HIJING
  primaries → real; Geant4 secondaries |id|==13 & IsPrimary==0 → still hadronic ✓).
- **② Parent-map barcode cutoff (`lim`=first bc>200000; map+BFS trace over `[0,lim)`): a SEPARATE,
  LEGITIMATE fix → KEEP.** Only fixes the HF-vs-prompt sub-split of already-real muons (barcode-
  duplication artifact, real-HF 5%→92%); does NOT touch real/hadronic. Aligns with "trace HF via
  Pythia truth only". Reverting would reintroduce the artifact.

**Plan (reviewer-gated):**
1. **Investigate** what the overlay index≥lim "untraced real" muons are (confirm HIJING per user's
   conditional label). Small diagnostic macro on the overlay NTUP mirroring the exact selection.
2. **d0_discrimination.C + dpop_dist.C**: delete the index-gate line; KEEP the parent-map cutoff;
   add a 5th overlay-only category **"real (untraced HIJING)"** (index<0 || index≥lim real muons;
   empty for signal since lim=ntr). → `/review-analysis-code` (quote this stage + §3).
3. **bkg_mc_provenance.C + fill_weighted_fullsim.C**: restore from `backup_pre_indexgate_20260707/`
   (100%-clean undo — no HF/prompt split there). Recompile, rerun fills, regen plots.
4. **Regenerate** d0/dpop plots → `/review-plot`. Verify signal identical (no-op), overlay real rises.
5. **Correct docs+memory**: doc E Progress Log/this stage, the sub-doc C note added during the
   scratch merge (currently describes the index gate as a fix — WRONG), memory
   `reference_muon_truth_provenance_nav`, and add the HIJING-real physics to the umbrella doc.

Data-file note (unrelated): `bkg_mc_provenance` `data_dos_dss_overlay.png` changed ~0.73× only because
its input `histograms_real_pairs_pp_2024_..._template_fit.root` was refreshed 2026-07-07 01:02 (backup
plot 2026-06-25) — a legit data update, not the classifier.

---
**(superseded — was CERTIFIED, now reopened for the index-gate revert) ✅ 2026-07-07.** Final `/review-plot` PASS (log
`review-plot-20260707-120911-upfront-d0-dpop-certify.md`, iter 2, 0 CRITICAL/0 WARNING; C1–C6 incl.
C5 provenance + C6 MC-weight). All 32 panels (both variables × 4-way primary + combined 3-way ×
signal + overlay) correct; both SUMMARYs match the run logs. Fixes this round: INDEX gate for "real"
(`muon_truth_index∈[0,lim)`; overlay real 91.7% HF / 0 unresolved; pp24 no-op), per-event MC weights
on BOTH modes (dpop mode-1 fixed; C6 enforcement added `b2e983d`), recolor red/blue/kGreen+2/magenta,
combined 3-way variant, descriptive sample labels/filenames (`fullsim_pp24`/`hijing_overlay`), Δp/p
cut-line label fixed. **Remaining:** (a) user's selection-tightening decision (Remaining Work). (b)✅ cosmetics resolved
(d0 subtitle moved to a centered top-margin title; dpop labels nudged in — both re-plotted &
verified). (c)✅ index gate applied to `bkg_mc_provenance_20260624` (2026-07-07 Progress Log): both
classifiers gated, signal EXACT no-op, overlay real→hadronic (prompt 96.57%→94.37%, had 3.05%→5.25%);
backup in `backup_pre_indexgate_20260707/`; a read-only agent is diffing backup-vs-new plots to
confirm no adverse effect on the data-vs-fullsim view. Plot dirs are data-area (not git).

---
**(superseded) ✅ 2026-07-07 — earlier PASS (before the index-gate + MC-weight round).**
Log `.claude/logs/review-plot-20260707-004523-upfront-bkg-d0-dpop-final.md` (iter 1, 0 CRITICAL/0
WARNING; provenance C5 + overlay-HF-recovery critical checks PASS; all numbers verified). Three
corrections that session, all applied to BOTH macros + regenerated + re-reviewed:
1. **OVERLAY provenance bug fixed.** Pre-fix overlay had real-HF ~5% + a huge fake "real-prompt
   residual" — an artifact of the HIJING-overlay **barcode duplication** (dups every event; Geant4/
   HIJING at barcode>200000). My standalone macros re-traced truth parents over the FULL merged truth
   record WITHOUT the generator-block cutoff → signal-HF muons' parents mis-mapped to colliding light
   hadrons. NOT a sample problem, NOT the mechanism failing: the cutoff (`PythiaTruthExtras.c:357`
   `pythia_only_barcode_cache` = first barcode>200000) lives only in the PAIR MC path; single-muon
   trees carry no MC truth origin, so a per-muon standalone macro must replicate it. Fix: build the
   barcode→index map + BFS trace only over `[0,lim)`, lim=first barcode>200000 → **overlay real-HF
   5%→90%** (matches signal 92%; displaced med 0.088 mm). No-op for signal. This is the case that
   motivated the new BLOCKING NTuple-Processing-Provenance rule (CLAUDE.md + `conventions/ntuple-
   provenance.md` + physics-results-review C5), committed `732b898`.
2. **Recolor** — 4 classes red / blue / kGreen+2 / magenta (real HF / real prompt / hadronic / fake).
3. **Combined 3-way variant added** (`combined/` subdir in each plot dir): real (HF+prompt) /
   hadronic / fake — user-requested supplementary purity view; 4-way stays primary. Composition of
   "real prompt" (the merged-in part): ~6–9% of real muons, dominated by light-meson resonances
   (ρ/ω/φ/η/η'→µµ) + quarkonium (J/ψ/ψ'/Υ) — all prompt (small |d0|) and downstream-separated in the
   minv fit as resonance peaks (veto / §3h templates), NOT part of the single-b S+G continuum.

**Tasks complete:** T1 WP doc, T2 binning, T3 |d0| (4-way + combined), T4 Δp/p (4-way + combined),
all `/review-plot` PASS. Plot dirs are data-area (not git). Awaiting the user's selection-tightening
decision (Remaining Work).

---

**✅ 2026-07-06 — (superseded by the reopen above for overlay) both plot sets PASSED `/review-plot`
(iter 2, 0 CRITICAL/0 WARNING).** Log `.claude/logs/review-plot-20260706-200246-upfront-bkg-d0-dpop.md`.
- **T1 (working points)** ✅ — `Analysis/docs/references/muon_working_points.md` (Tight adds χ²<8 +
  binned ρ′/(q/p)-sig; graded recommendation). Committed `80de5a4`.
- **T2 (single-μ pT binning)** ✅ — `single_mu_pt_coarse_bins={4,8,14,25,100}` in ParamsSet.h.
  Committed `cec2416`.
- **T3 (|d0| plots)** ✅ — `d0_discrimination_20260706/` (8 PNGs + SUMMARY.md). Classification bug
  (truth_vtx-capped ntr) found+fixed; BFS chain-walk nav. Physics: signal real-prompt small-|d0|,
  real-HF displaced; overlay hadronic SMALLEST |d0| (0.021 mm) → hypothesis confirmed.
- **T4 (Δp/p plots)** ✅ — `dpop_distribution_20260706/` (8 PNGs + SUMMARY.md). Δp/p real symmetric,
  hadronic broad+positive; ~50–60% hadronic survives Δp/p<0.12 at all pT (mild boost rise) → partial
  discriminant, motivates |d0|.
- **Doc factorization** ✅ — this is sub-doc E; umbrella + A–D committed `522484e`.
Plot dirs are data-area (not git). This doc + review log are the git-tracked record.
**Clear Latest Stage — no active thread; awaiting the user's selection-tightening decision above.**
