# Muon Working-Point (WP) Registry — where the Medium/Tight choice lives

**Status:** LIVING reference. **Default WP: TIGHT** — IMPLEMENTED & RERUN 2026-07-08 (changed from
Medium; code committed `11b0748`, `/review-analysis-code` PASS; tight trig-eff turn-on + tight reco-eff
+ tight crossx/R_AA regenerated and `/review-plot` PASS; see `docs/tracking/tight_wp_default_change.md`).
**One open item:** pp-tight reco-eff has no source → INTERIM reuse of Medium Fig.31 (labeled); real
options = keep interim OR peripheral-PbPb-tight proxy (user decision; §4). **Purpose:** single place listing every site where the
muon quality WP is chosen/applied/labeled, so the **WP systematic uncertainty** can be evaluated by
switching Tight→Medium and quantifying the shift in the final result. **Maintenance:** update this
table whenever a new WP-sensitive place is added as the analysis is completed.

## WP definition
`muon_quality` bitmask: **MEDIUM = `& (1|8|32|256)`** (combined|Medium|IDCuts|MuonCuts);
**TIGHT = `& (1|16|32|256)`** (bit 8→16; Tight ⊂ Medium since `getQuality()` Tight=0<Medium=1). Only
this bit changes for the systematic. Full definition: `docs/references/muon_working_points.md`.
There is **NO single source of truth** — the choice is duplicated across the sites below (a future
refactor could centralize it in `ParamsSet.h`).

## Systematic procedure
Nominal = Tight. To evaluate the WP systematic: set every "knob" site below to Medium, rerun the
affected chain (see `signal_selection_change_impact.md` + the blast-radius notes here), and take the
difference in the final observable (crossx / R_AA) as the WP uncertainty term. **Result-affecting**
sites (marked ⚠) move the number; **cosmetic** sites (labels) do not but must be kept consistent.

## Registry (by pipeline stage)

### 0. Skim — NOT a knob (both bits already stored)
- `SkimCode/source/HFtrigValidation/src/TrigRates.cxx:1196-1197` — packs BOTH `+=8` (Medium) and
  `+=16` (Tight) into `muon_quality`. ⇒ **switching WP never needs a re-skim.** (`MyUtils.h:10-11`
  `HI_TIGHT`/`HI_LOOSE` are track quality for centrality — NOT the muon WP; exclude.)

### 1. NTuple processing — MC ⚠ (both flags written; select downstream — no reprocess)
- `NTupleProcessingCode/PythiaFullSimExtras.c:34-37` `PassMuonMediumCuts` (`&1&8&32&256`);
  `:136-137` `pass_medium`/`pass_tight`(`&16`); `:240-241,283-286` `pair_pass_medium`/`pair_pass_tight`
  (+`_and_resonance`). Both flags written unfiltered.
- `NTupleProcessingCode/PowhegFullSimExtras.c:29-32, 154-155, 206-207, 244-245` — same.
- `PythiaFullSimOverlayExtras.c` (PbPb overlay) inherits `PythiaFullSimExtras` masks.
- Struct fields: `MuonObjectsParamsAndHelpers/Muon.h:18,39`, `MuonPairMC.h:53-55`, `MuonPairReco.h:29`.

### 2. NTuple processing — DATA ⚠ (WP is a HARD CUT here; default was Medium)
- `NTupleProcessingCode/DimuonDataAlgCoreT.c:578-590` `PassCuts_DataCore(bool requireTight)`:
  `&1`; `if(requireTight) &16 else &8`; `&32`; `&256`. Called `:673`.
- `:396` `tight_suffix = requireTight?"_tight":""`; `:675` `pair_pass_tight` tagged (`&16`) regardless.
- `DimuonDataAlgCoreT.h:340` `bool requireTight = false;` ← **the DATA WP default**.
- Driver knob (commented): `NTupleProcessingCode/run_pp_24.sh:21`.
- **Two routes to Tight (no re-skim either way):** (a) reprocess data with `requireTight=true`
  (regenerates `muon_pairs_*` trees, `_tight` suffix), OR (b) **preferred, no reprocess:** wire the
  dead RDF `isTight` (see §3) to `Filter("pair_pass_tight")` — the flag is already stored (`:675`).
- Legacy `original_no_template_class/DimuonDataAnalysisBaseClass.*` — same logic, inactive.

### 3. RDF hist-filling ⚠
- **DATA:** `RDFBasedHistFillingData.h:176` `bool isTight=false;` — **declared but UNUSED** (`.cxx:32`
  help-print only; no `Filter`). `RDFBasedHistFillingPP.cxx:545` / `RDFBasedHistFillingPbPb.cxx`
  inherit the ntuple WP (no re-apply). ⇒ to select Tight without reprocessing, wire `isTight` →
  `Filter("pair_pass_tight")` here.
- **MC (already both variants):** `RDFBasedHistFillingPythiaFullsim.cxx`, `...Overlay.cxx` (PbPb),
  `...PowhegFullsim.cxx`, `...PowhegFullsimSingleMuon.cxx` build both `_pass_medium` & `_pass_tight`
  hist variants — just select the tight output.

### 4. Reco-efficiency ⚠ (Tight input EXISTS; active placeholder built from MEDIUM)
- Active placeholder: `plotting_codes/reco_effcy/build_run2_reco_eff_placeholder.C:85` reads
  `EffFiles/MuonRecoEffcyRun2MC_medium.root` (Tight analog `MuonRecoEffcyRun2MC_tight.root` EXISTS but
  is not used); output keys/labels `:107,121-163` say "medium".
- Consumed in crossx: `RDFBasedHistFillingData.cxx:656` loads `run2_reco_eff_placeholder.root`;
  `:698-699,715` **hardcode medium keys** (`gr_reco_eff_medium_pp_*`, `tf1_reco_eff_medium_pbpb_*`).
  ⇒ **Tight default REQUIRES rebuilding the placeholder from `_tight.root` + switching these keys** —
  else a Tight-selected spectrum gets a mismatched Medium reco-eff correction.
- PbPb: `EfficiencyCorrs/NoOverlayMC.C/.h:126,720-726,870` (default MEDIUM).
- Reco-eff plotters (have a `tight_WP` flag, default false=medium; drivers default medium):
  `PythiaFullsimRecoEffPlotter.cxx:28-29,110,131,146`, `PowhegFullsimRecoEffPlotter.cxx:30-32,175,...`,
  `PowhegFullsimDetRespPlotterSingleMuon.cxx:24-26,86,...`.
- **HARDCODED medium (needs a WP config var):** `plotting_codes/reco_effcy/plot_single_muon_reco_effcy.cxx`
  (`Filter("pass_medium")` `:131`, y-title `:169,229,317`) and `..._r17618_vs_r17662.cxx`.
- Legacy JPsi path: `EfficiencyCorrs/Bins.h:480-522` (`Reco_Medium/Tight_JPsi.root`; files absent).

### 5. Trigger-efficiency ⚠ (framework supports both, default MEDIUM; active pipeline = implicit medium)
- `EfficiencyCorrs/Bins.h:146-162` `enum QualityCut{MEDIUM=0,TIGHT=1}` + `PassQualityCut` (bit8/bit16)
  + labels; `:640` default `m_quality_cut = MEDIUM`; `:724-748` `"_Tight"` name suffix.
- `EfficiencyCorrs/TrigEff.C:96,282-286,541,548` (default MEDIUM); `TrigEffMu4NoL1.C:127,140,1017`
  (default MEDIUM, `l_qual=8`).
- `EfficiencyCorrs/TrigAndRecoEff.C:198-203` — **no direct Tight trigger-eff**; Tight = Medium ×
  (Tight/Medium ratio). ⇒ a genuine Tight default needs the Tight trig-eff derived/validated.
- Active PbPb pipeline `plotting_codes/trig_effcy/TrigEffPlotterPbPb.cxx:561` applies **no WP filter**
  (implicit medium via input trees); label hardcodes `"Medium #mu"`.

### 6. Plotting labels — cosmetic (must stay consistent; will go stale)
`plot_single_b_crossx_pp.cxx:28,33,45,50`; `plot_single_b_crossx_pbpb.cxx:14,30`;
`plot_single_muon_reco_effcy.cxx:169,229,317`; `plot_reco_distr_singleb_vs_op_pp24.C:16,88,89,120,121`;
`TrigEffPlotterPbPb.cxx:561`; `build_run2_reco_eff_placeholder.C:151,163`;
`dphi_plots/dphi_mc_data_compr_with_tight.c` (legacy dphi study).

### 7. Study macros (data area, not git-tracked)
- `dimuon_data/plots/template_fitting/d0_discrimination_20260706/code/d0_discrimination.C`,
  `dpop_distribution_20260706/code/dpop_dist.C` — `QBITS` (currently Medium; need a `useTight` config
  var, default Tight — per memory `feedback_plots_wp_config_var`).
- `bkg_mc_provenance_20260624/code/{bkg_mc_provenance.C,fill_weighted_fullsim.C}` — `useTight` already present.

## Blast radius of a Medium↔Tight switch (see also `signal_selection_change_impact.md`)
That doc's "ntuple UNCHANGED" assumption is FALSE for the WP (WP is upstream for DATA). Net
result-affecting chain: DATA WP selection (§2/§3) → data crossx (pp + PbPb 23/24/25) → crossx / R_AA /
MC–data / stage plots; **reco-eff placeholder must be rebuilt from the Tight file (§4)** and
**trig-eff must be Tight (§5)** — both move the number (unlike a ΔR change where the placeholder
cancels). Truth acceptance is WP-independent (unchanged). Skim unchanged. MC RDF already has both.

---
*Referenced from the analysis roadmap/status for systematic-uncertainty evaluation. Keep in sync with
`docs/tracking/tight_wp_default_change.md` (the implementation) and memory
`feedback_plots_wp_config_var` (all plot sets must expose a Medium/Tight config var, default Tight).*
