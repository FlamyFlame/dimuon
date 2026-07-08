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

## Latest Stage
**✅ Stage 1 (Discovery) COMPLETE — registry produced.** **⏸ Stage 2: PRESENTING THE PLAN TO THE USER**
before executing the result-affecting parts (reco-eff-Tight placeholder rebuild, Tight trig-eff, data
crossx reruns for pp + PbPb 23/24/25). Awaiting user direction on execution scope/sequencing. The
tight d0/dpop study-plot set (add `useTight`, default Tight) is the low-risk first execution step.
