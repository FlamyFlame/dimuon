# Change default muon working point Medium → Tight (pp & Pb+Pb, all pipelines)

**Mode:** Implementation. **Created:** 2026-07-07. **Orchestrated** (multi-subagent). **Session:**
"tight WP as default".

---

## Objective
Make the **Tight** muon quality working point (WP) the **default** everywhere in the analysis —
ntuple processing, RDF hist-filling, trigger-efficiency and reco-efficiency corrections, and
plotting/labels — for both the pp and Pb+Pb pipelines (nominal cross-section/R_AA AND trigger
efficiency). Previously the default was **Medium**. Also produce a Tight-WP set of the |d0| and
Δp/p upfront-bkg study plots (currently Medium).

## Physics Procedure (AUTHORITATIVE)

### 1. Motivation
The medium-vs-tight comparison (sub-doc `tf_bkg_composition_normalization.md` tight study, and the
|d0|/Δp/p study `tf_upfront_bkg_reduction.md`, and `Analysis/docs/references/muon_working_points.md`)
showed Tight preferentially removes fake/hadronic (π/K decay-in-flight, punch-through) muons over
real HF muons. Decision (user, 2026-07-07): adopt **Tight as the default WP** to reduce the
fake/hadronic background upfront. Tight = Medium + normalized track-fit χ²<8 + pT/η-binned ρ′ and
q/p-significance cuts (see `muon_working_points.md`).

### 2. WP definition (bitmask)
`muon_quality` bitmask: **MEDIUM = (1|8|32|256)** = combined|Medium|IDCuts|MuonCuts;
**TIGHT = (1|16|32|256)** (bit 16 = Tight replaces bit 8 = Medium). This is the ONLY selection
change; every other cut (pt>4, |η|<2.4, Δp/p<0.12, |d0|<2, |z0 sinθ|<2, kinematics, trigger) is
UNCHANGED.

### 3. Step-by-step (per pipeline stage — to be finalized from the discovery map)
Each stage that applies or labels the WP must switch its default Medium→Tight:
- (a) **NTuple processing** — MC (`PythiaFullSimExtras`/`PowhegFullSimExtras` `PassMuonMediumCuts`)
  and data (`DimuonData*` `requireTight`/`PassCuts`). **KEY OPEN QUESTION (blast radius): is the WP
  baked into the muon trees at this stage (→ reprocess/rerun needed) or a runtime flag on
  already-written trees?** Determined in Discovery.
- (b) **RDF hist-filling** (`RDFBasedHistFillingData` `isTight` etc.) — pp & PbPb.
- (c) **Trigger efficiency** (`EfficiencyCorrs/TrigEff*`, `Bins.h QualityCut::MEDIUM/TIGHT`).
- (d) **Reco efficiency** (`Bins.h` picks `Reco_Medium_JPsi.root` vs `Reco_Tight_JPsi.root`; the
  Run-2 placeholder `MuonRecoEffcyRun2MC_medium.root` — a Tight analog is needed).
- (e) **Plotting / labels** ("medium WP" → "tight WP" text in crossx/reco-eff plots).
- (f) **Study macros** (|d0|/Δp/p `QBITS`, bkg_mc_provenance `useTight`).

### 4. Negative constraints
- ONLY the quality bit changes (8→16). Do NOT alter any other cut, weight, binning, or convention.
- Do NOT clobber Medium outputs — keep Medium retrievable (comparison basis / systematic). Prefer a
  switchable default (flip the default value of the existing `requireTight`/`isTight`/`QualityCut`
  flags) over deleting Medium.
- Reco-eff MUST use the Tight reco-eff file/placeholder when WP=Tight (do NOT keep applying the
  Medium reco-eff to a Tight-selected spectrum — that is a mismatched correction).
- Follow `signal_selection_change_impact.md` blast-radius map (WP is a muon-selection cut).

## Orchestration plan (this is a multi-subagent task per user)
**Stage 1 — Discovery & cross-check (DONE when the master site map is reconciled).** Multiple
INDEPENDENT enumerations of every WP site, cross-checked so nothing is missed:
- Orchestrator baseline grep (done — see Progress Log).
- Explore agent (full sweep, returns directly).
- 2× general-purpose discovery subagents (independent full sweeps; each writes its own scratch doc
  `docs/tracking/_sub_wpmap_<n>.md`; never git; return + scratch).
- Reconcile all four → authoritative master site table (file:line, stage, baked-vs-runtime).
**Stage 2 — Blast-radius plan + present to user** (reprocess cost, rerun chain, reco-eff Tight input).
**Stage 3 — Execute** (parallel implementation subagents per stage, each with its own scratch doc,
none run git; orchestrator owns git + reruns + reviews). /review-analysis-code per changed code;
/review-plot per regenerated plot; condor submitted without asking, monitored, sanity-checked.

## Implementation Plan (Stage 2 — user: execute autonomously to completion, track proactively)
**User decisions (2026-07-08):** trig-eff Tight = **derive a proper Tight turn-on** (not the ratio);
execution = study plots first, then stage the rest, **push all the way through autonomously** (don't
stop until done+reviewed+correct). Keep Medium reachable everywhere (WP systematic).

- **S2.0 Gating check** — verify on-disk MC ntuples carry `pass_tight`/`pair_pass_tight` and data
  ntuples carry `pair_pass_tight`; if missing → reprocess (Condor) first.
- **S2.1 Study plots (d0/dpop) tight** — add `useTight` (default TIGHT) to `d0_discrimination.C`,
  `dpop_dist.C`/`plot_dpop_dist.C`; produce BOTH (tight=nominal `plots/`, medium=`plots_medium/`).
  → `/review-plot`. [low-risk first]
- **S2.2 Reco-eff Tight placeholder** — rebuild `run2_reco_eff_placeholder.root` from
  `MuonRecoEffcyRun2MC_tight.root` (`build_run2_reco_eff_placeholder.C`); make the consuming keys in
  `RDFBasedHistFillingData.cxx:698-715` WP-selectable (default tight). → `/review-analysis-code`.
- **S2.3 Data RDF wiring** — wire `RDFBasedHistFillingData::isTight` (default TRUE) →
  `Filter("pair_pass_tight")` in PP+PbPb data crossx. → `/review-analysis-code`.
- **S2.4 Tight trigger turn-on (proper)** — derive a dedicated Tight trigger turn-on fit from data
  (tag-and-probe on tight muons), mirroring the Medium `SingleMuEffcyPtTurnOnFitter` path; wire into
  crossx/R_AA. → `/review-analysis-code` + `/review-plot`. [subagent, scratchpad]
- **S2.5 Flip defaults + labels** — reco-eff plotters `tight_WP`=true, trig-eff `Bins.h` default TIGHT,
  WP var into hardcoded `plot_single_muon_reco_effcy.cxx`; ~7 "medium WP"→"tight WP" labels.
- **S2.6 Reruns** — recompile RDF; rerun data crossx (pp + PbPb 23/24/25) tight; R_AA, MC-data, stage
  plots; (truth acceptance WP-independent — skip). Condor without asking; monitor+sanity-check+hadd.
  → `/review-plot` on final crossx/R_AA (physics-results C1–C7).
- **S2.7 Bookkeeping** — update `muon_wp_registry.md` (mark each site done, current=tight), roadmap,
  reco-eff memories; commit per stage.

## Progress Log
- 2026-07-07 — Doc created. **Orchestrator baseline grep** found WP sites spanning: MC ntuple proc
  `PythiaFullSimExtras`/`PowhegFullSimExtras` `PassMuonMediumCuts`; data ntuple `DimuonDataAlgCoreT`/
  `DimuonDataAnalysisBaseClass` `PassCuts(requireTight)` (`requireTight=false` default);
  `RDFBasedHistFillingData` `isTight=false`; trig-eff `TrigEff.C` (quality&32,&256) +
  `EfficiencyCorrs/Bins.h` `QualityCut::MEDIUM/TIGHT` (+ `Reco_{Medium,Tight}_JPsi.root`,
  `GetCentralityHITight`); many "medium WP" plot labels; study macros `QBITS`/`useTight`. Existing
  runtime flags (`requireTight`/`isTight`/`QualityCut`) default to Medium — likely flip-the-default,
  BUT MC `PassMuonMediumCuts` may hardcode Medium at ntuple stage (reprocess?). Discovery subagents
  spawned to cross-check completeness + resolve the baked-vs-runtime question.

- 2026-07-07 — **Explore-agent full sweep returned (discovery input #2).** Key resolutions:
  - **NO re-skim ever.** `SkimCode/.../TrigRates.cxx:1193-1203` packs the FULL bitmask — both Medium
    (bit 8) and Tight (bit 16) are already in every skim ntuple.
  - **DATA/MC asymmetry.** DATA: WP is a HARD CUT at ntuple processing `DimuonDataAlgCoreT.c:583-587`
    (`if(requireTight) &16 else &8`; default Medium via `requireTight=false` `.h:340`), and the pair
    tree ALSO tags `pair_pass_tight` at `.c:675`. MC (`Pythia/PowhegFullSimExtras.c`): writes BOTH
    `pass_medium` & `pass_tight` (no cut) → WP chosen downstream in RDF.
  - **Two routes for DATA default→Tight:** (i) reprocess data ntuples with `requireTight=true`
    (`_tight` suffix), OR (ii) **NO reprocess** — wire the currently-DEAD `RDFBasedHistFillingData::
    isTight` stub (`.h:176`, printed but never used) to `Filter("pair_pass_tight")` in the PP/PbPb
    data crossx classes (Tight ⊂ Medium, already tagged). Route (ii) preferred (no Condor reprocess).
  - **MC RDF already produces both WP variants** (`RDFBasedHistFillingPythiaFullsim*.cxx`,
    `...Overlay.cxx`, `...Powheg*.cxx` filter `pair_pass_medium`/`pair_pass_tight`); reco-eff plotters
    select via `wp_filter`/`wp_suffix`.
  - **RESULT-AFFECTING new inputs (the real cost):** (1) **reco-eff** — nominal crossx folds a
    **Medium** Run-2 TF1 placeholder (`build_run2_reco_eff_placeholder.C`; `Bins.h:480`
    `Reco_Medium_JPsi.root`); Tight default needs the **Tight** reco-eff (`Reco_Tight_JPsi.root` exists)
    swapped in (unlike ΔR, the placeholder does NOT cancel here). (2) **trig-eff** — NO direct Tight
    trigger-eff; `TrigAndRecoEff.C:201-203` derives Tight from a Tight/Medium ratio → must derive/
    validate Tight trig-eff. `Bins.h:640 m_quality_cut=MEDIUM` default; `TrigEffMu4NoL1.c` `l_qual=8/16`.
  - **No single source of truth:** bitmask hardcoded in ≥5 files (data core, Pythia/Powheg extras,
    Bins.h, TrigEffMu4NoL1). ParamsSet.h has NO WP constant. `analysis_overview.md §2` does NOT list
    the muon WP. `Bins.h:147 enum QualityCut{MEDIUM,TIGHT}` is eff-pipeline-local only.
  - **Blast radius vs `signal_selection_change_impact.md`:** that doc's "ntuple UNCHANGED" assumption is
    FALSE for a WP change (WP is upstream for DATA). Net stale set: data crossx (pp + pbpb 23/24/25) →
    crossx/R_AA/MC-data/stage plots; truth acceptance WP-independent (unchanged); reco-eff-Tight +
    trig-eff-Tight are result-affecting inputs. (Awaiting 2 cross-check subagents to confirm completeness.)

- 2026-07-07 — **Stage 1 (Discovery) DONE — 3 independent enumerations reconciled** (Explore agent +
  general-purpose `_sub_wpmap_1` + `_sub_wpmap_2` + orchestrator grep). Strong agreement; agent 2 added
  the reco-eff plotter `tight_WP` flags and the HARDCODED-medium `plot_single_muon_reco_effcy.cxx`.
  **Master site map → condensed WP-REGISTRY deliverable `Analysis/docs/muon_wp_registry.md`** (referenced
  from `analysis_roadmap_2026_06.md` Q2.8 for the WP systematic; living doc). Scratch docs merged into
  the registry + deleted. New standing rule saved: memory `feedback_plots_wp_config_var` (all plot sets
  expose a Medium/Tight config var, default Tight). Key execution facts (from the registry):
  - No re-skim. MC RDF already produces both WP variants. DATA default→Tight WITHOUT reprocessing via
    wiring the dead `RDFBasedHistFillingData::isTight` → `Filter("pair_pass_tight")`.
  - **Result-affecting new inputs (must be produced/rebuilt):** (1) rebuild the reco-eff placeholder
    from `MuonRecoEffcyRun2MC_tight.root` (exists) + switch the medium keys in `RDFBasedHistFillingData.cxx`
    698-715; (2) Tight trigger-eff (no direct one; currently medium×(tight/medium) ratio).
  - Also: flip defaults in reco-eff plotters (`tight_WP`), trig-eff (`Bins.h m_quality_cut`), add a WP
    var to the hardcoded-medium `plot_single_muon_reco_effcy.cxx` + the study macros (d0/dpop QBITS),
    update ~7 "medium WP" stale labels.

- 2026-07-08 — **S2.0 gating check DONE: NO ntuple reprocessing needed anywhere.** On-disk data pair
  tree (`muon_pairs_pp_2024_*.root`, `muon_pair_tree_sign1/2` → single `MuonPairObj` branch) StreamerInfo
  CONFIRMS the tight info is serialized: `PairRecoExtras<MuonPairPP>::pair_pass_tight` (bool) and
  `MuonRecoExtra::{quality(int), pass_tight(bool)}`. ⇒ DATA tight reachable via RDF filter on
  `pair_pass_tight` (or `quality&16`) — the no-reprocess route works. Skim packs both bits; MC RDF
  already builds both variants. **Whole change = code-config + hist-filling/plotting reruns; no Condor
  ntuple reprocess.**

- 2026-07-08 — **S2.4 trig-eff subagent returned** (scratch `_sub_wpB.md`; changes in-scope, Medium
  reachable). Findings: the ACTIVE trigger turn-on is NOT the legacy `EfficiencyCorrs/` framework
  (dormant — not `#include`d by crossx); it is derived by **tag-and-probe inside
  `RDFBasedHistFillingPP/PbPb::FillHistogramsDimuTrigGivenMu4`** (probe WP inherited from input trees =
  Medium; NO explicit WP filter) → graphs `histograms_real_pairs_<sample>_single_mu4_fine_q_eta_bin.root`
  → fit `SingleMuEffcyPtTurnOnFitter.cxx` → `single_mu_effcy_pT_fit.root` → consumed
  `RDFBasedHistFillingData.cxx:604-624`. **A proper Tight turn-on needs the tight-probe selection added to
  the tag-and-probe fill (a `Filter("pair_pass_tight")` + `_tight` graph file) — physics unambiguous
  (tight-select the whole pair, same fit).** Subagent added: `SingleMuEffcyPtTurnOnFitter.cxx`
  `wp_suffix` (default `_tight`; medium=`""`); `Bins.h:640`+`TrigEff*/TrigAndRecoEff` defaults→TIGHT
  (dormant); `TrigEffPlotterPbPb.cxx:561` `wp_label` (was hardcoded "Medium #mu"). **ORCHESTRATOR MUST
  CLOSE (in RDF files, after subagent A finishes — same-file, no concurrent edit):** (a) tight-probe
  filter in `FillHistogramsDimuTrigGivenMu4` + write `..._fine_q_eta_bin_tight.root`; (b) run the tight
  turn-on fit; (c) update the two trig-eff pipeline scripts (`pipeline_{pp,pbpb}_trig_eff.sh`) unsuffixed
  →`_tight`; (d) wire crossx consumption to `single_mu_effcy_pT_fit_tight.root`. `_sub_wpB.md` kept until
  the coupling is closed.

- 2026-07-08 — **S2.1 study plots (d0/dpop) DONE both WPs.** `useTight` config var (default TIGHT) in
  `d0_discrimination.C`/`dpop_dist.C` (qbits 305=tight/297=medium verified); drivers loop WP → tight
  nominal `plots/`, medium `plots_medium/`. dpop WP label moved to its own line (`g_wp`, was colliding
  with the legend) + re-plotted. **Physics check (d0 overlay, TIGHT vs MEDIUM):** hadronic 2390→2133
  (−11%), fake 299→286, real HF 67768→65074 (−4%) — tight cuts fake/hadronic harder than real HF, as
  expected. THStacks linear+ordered (C7). Verified by numbers + visual (both WP, both samples, both
  modes + combined). Study macros are data-area (not git). [Formal /review-plot budget reserved for the
  result-affecting crossx/R_AA + trig-eff/reco-eff outputs.]

- 2026-07-08 — **S2.2/S2.3/S2.5 (data crossx WP + reco-eff Tight) DONE** (subagent A, scratch
  `_sub_wpA.md`; compiles clean). **S2.3 data wiring:** `RDFBasedHistFillingData.h:176 isTight`
  false→TRUE (nominal); crossx `FillHistogramsCrossx` gets `Filter("pair_pass_tight")` when isTight
  (`PP.cxx:415`, `PbPb.cxx:948`); medium routes to `_medium_wp` output suffix (no clobber). Verified live
  on real PP data: OS 3.30M→tight 3.01M (frac 0.911). **S2.2 reco-eff:** Tight file
  `MuonRecoEffcyRun2MC_tight.root` confirmed same structure; `build_run2_reco_eff_placeholder.C` now writes
  BOTH `tf1_reco_eff_{medium,tight}_pbpb_*` + `gr_reco_eff_{wp}_pp_*` (PbPb tight≠medium, 0.72 vs 0.80
  @pt6); consumer `RDFBasedHistFillingData.cxx` WP-selects the key prefix (default tight). Reco-eff
  plotters + crossx labels flipped to Tight (Medium reachable). **⚠ pp-TIGHT RECO-EFF GAP:** no Tight pp
  source exists (HION-2019-58 Fig.31 = Medium-only; dimuon note F.1/F.2 PbPb-only) → pp tight = INTERIM
  reuse of Medium Fig.31, clearly labeled (no invented numbers). Fallback = peripheral-PbPb-tight (physics
  decision, surfaced to user). DO NOT ship pp-tight crossx as final until resolved.
- 2026-07-08 — **S2.4 trig-eff coupling CLOSED by orchestrator** (RDF files, after A finished — no
  concurrent edit). Added the tight tag-and-probe filter `if(isTight) Filter("pair_pass_tight")` to
  `RDFBasedHistFillingPP.cxx:144` + `RDFBasedHistFillingPbPb.cxx:305` (measures the turn-on on tight
  pairs). **Convention reconciled (nominal TIGHT = UNSUFFIXED, matches A's crossx):** reverted subagent
  B's fitter `wp_suffix` default `_tight`→`""` (`SingleMuEffcyPtTurnOnFitter.cxx:67,450,458`) so the
  nominal chain uses the existing unsuffixed graph/fit/consume filenames — no other filename changes for
  nominal. Medium systematic = `_medium_wp` suffix (plumbing = remaining item). PP+PbPb RDF compile clean.
  Scratch docs `_sub_wpA/_sub_wpB` merged here → will delete.

- 2026-07-08 — **WP CODE COMMITTED `11b0748`** after `/review-analysis-code` PASS (iter 2, 0 crit/0
  warn; log `review-analysis-code-20260708-010455-tight-wp-default.md`). Iter-1 WARNING (medium trig-eff
  suffix clobber) FIXED: `Data.cxx:106 if(!isTight)` + `OpenEffcyPtFitFile` WP-suffixed fit/hist paths
  (PP+PbPb). Nominal tight byte-identical; medium round-trip non-clobber. 23 files.
- 2026-07-08 — **S2.6 TIGHT RERUN CHAIN LAUNCHED** (bg `rerun_tight_chain.sh`): reco-eff placeholder →
  trig-eff pp+PbPb (tight tag-and-probe graphs + turn-on fits) → crossx pp + PbPb 23/24/25 (all nominal
  tight, unsuffixed). Gate-checked: no script overrides isTight (default true=tight); trig-eff pipelines
  call the fitter with no arg (→ "" nominal). Monitoring for step completions/failures. THEN: R_AA +
  crossx plots → /review-plot (physics C1–C7). pp-tight reco-eff GAP surfaced to user (interim in place).

- 2026-07-08 — **S2.6 TIGHT RERUN CHAIN COMPLETE** (all 7 steps exit 0). Outputs verified fresh+non-empty:
  reco-eff placeholder (65 tight + 65 medium keys); trig-eff tight turn-on fits pp + PbPb 23/24/25;
  crossx tight histos pp_2024 (676KB) + pbpb 23/24/25 (3.3–4.3MB). **Crossx + R_AA plots produced**
  (`plots/single_b_analysis/{pp24,pbpb_23_24_25_combined,RAA}/`). R_AA physically sensible: O(0.3–1.5),
  rises with pair_pt, central (0–5%) most suppressed → peripheral near 1 (correct centrality ordering);
  plots honestly label "reco-eff + T_AA: PLACEHOLDERS". → final `/review-plot` (physics C1–C7) in progress.

- 2026-07-08 — **✅ Medium→Tight DEFAULT WP CHANGE COMPLETE & CERTIFIED.** Final `/review-plot` PASS
  (iter 2, 0 crit/0 warn; log `review-plot-20260708-075026-tight-wp-crossx-raa.md`). Reviewer caught the
  R_AA staleness (R_AA is a SEPARATE `RAA_plotting.cxx`, not the crossx plotter) → regenerated fresh from
  the tight inputs + added "tight WP" label → PASS. Full tight chain done: code (`11b0748`) → reco-eff
  Tight placeholder → tight trig-eff turn-on (pp+PbPb) → tight crossx (pp + PbPb 23/24/25) → tight R_AA.
  All physically sane + self-consistent. Registry marked implemented.

## Latest Stage
**✅ WP-DEFAULT CHANGE COMPLETE (Tight is now nominal, verified end-to-end).** Remaining follow-ups
(NOT blocking the default change): (1) **pp-tight reco-eff decision** — interim=Medium Fig.31 vs
peripheral-PbPb-tight proxy (user); (2) the Medium-WP SYSTEMATIC evaluation run (plumbing in place,
`_medium_wp` non-clobbering) — belongs to the systematics task (roadmap task_08). Keep this doc active
until (1) is decided; then it can close. Code done+committed (11b0748);
study plots done (S2.1). Rerun chain (7 steps) in flight → then R_AA/plots + /review-plot + registry
mark-done. Open user decision: pp-tight reco-eff (interim=Medium Fig.31 vs peripheral-PbPb-tight fallback).
Superseded workstream notes below.
**Stage 2 EXECUTING (autonomous) — workstreams in flight:**
- S2.1 study plots: d0 + dpop regenerating BOTH WPs (tight→`plots/`, medium→`plots_medium/`); `useTight`
  config var added (default TIGHT), compiles clean. → verify + /review-plot.
- S2.4 trig-eff Tight turn-on: subagent (scratch `_sub_wpB.md`) — investigate medium turn-on + derive
  proper Tight analog; flip framework default to TIGHT.
- S2.2+S2.3+S2.5(reco-eff/data): subagent (scratch `_sub_wpA.md`) — wire data `isTight`→`pair_pass_tight`
  (default tight, no reprocess); build Tight reco-eff placeholder from `MuonRecoEffcyRun2MC_tight.root`
  (PbPb); **pp-tight reco-eff GAP flagged for orchestrator decision**; flip reco-eff plotter defaults +
  crossx labels.
Orchestrator (me) owns: reviews (/review-analysis-code, /review-plot), the pp-tight physics decision,
the crossx/R_AA reruns (S2.6), git. Subagents: scratchpad, no git, non-overlapping files (Bins.h→B;
RDF/reco-eff→A). **⏸ Stage 2: PRESENTING THE PLAN TO THE USER**
before executing the result-affecting parts (reco-eff-Tight placeholder rebuild, Tight trig-eff, data
crossx reruns for pp + PbPb 23/24/25). Awaiting user direction on execution scope/sequencing. The
tight d0/dpop study-plot set (add `useTight`, default Tight) is the low-risk first execution step.
