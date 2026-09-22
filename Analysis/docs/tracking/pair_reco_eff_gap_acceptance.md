# Gap-cut acceptance folded into the pair reconstruction efficiency

Mode: **IMPLEMENTATION**. Opened 2026-09-17.

## Objective

User decision (2026-09-17): the **detector-gap cuts** — the single-muon q·η fiducial windows
`ParamsSet::single_mu_fiducial_gap_cuts` applied to both muons and the pair-level
`|η^pair| < ParamsSet::pair_eta_fiducial_max = 2.2` — are **acceptance** cuts, imposed because the
detector gaps make the muon acceptance / reconstruction efficiency low and rapidly varying
there. They are NOT signal cuts (the mass window and the pair-pT threshold, which exist to
improve signal/background). Their cost is therefore to be carried by the **pair-level
reconstruction efficiency** ε_reco, not by the definition of the truth-level signal region and
not by a separate acceptance factor ε_acc. Implement that: the truth denominator of ε_reco drops
the gap cuts; the reco numerator keeps them (it must stay identical to the data selection).

Edge cases the user named explicitly: a pair with truth |η^pair| < 2.2 whose reconstructed
|η^pair| ≥ 2.2 is killed by the gap cut in data, so in MC it must sit in the denominator and NOT
in the numerator — it is an inefficiency.

This supersedes the "ε_acc, not yet built" item of `muon_gap_cuts_acceptance.md` (F12/F17): ε_acc
is no longer a separate factor.

## Physics Procedure

### 1. Motivation

The efficiency-corrected yield estimates the number of truth signal pairs in a truth fiducial
region. Every reconstruction-level requirement that a truth pair in that region can fail is an
inefficiency and belongs in ε_reco. Until now the gap cuts were applied on BOTH legs (truth
q·η / truth η^pair in the denominator, reco in the numerator), making ε_reco a *fiducial*
efficiency and the cross-section a *gap-fiducial* cross-section that still needed an unbuilt
truth-level factor ε_acc (`muon_gap_cuts_acceptance.md` F12: 0.88, stale). Folding the gap cuts
into ε_reco removes that factor, makes the truth fiducial region the physically simple one
(muon p_T > 4.5 GeV, |η| < 2.4; pair mass window, pair p_T > 9 GeV; |η^pair| < 2.2 as the
measured range — see §2), and makes ε_reco carry the gap loss **differentially** in (pair p_T, η^pair, ΔR) instead of as one global number.

### 2. Top-level equation

Unchanged in form (`analysis_overview.md` §4b):

    dσ/dX = (1 / L_int) · Σ_{reco pairs passing the DATA signal selection} w_reco · w_trig,
    w_reco = 1 / ε_reco(pair p_T, η^pair, ΔR)   evaluated at RECO kinematics (pre-unfolding)

**DATA signal selection** (unchanged; `RDFBasedHistFillingPP.cxx` `signal_cuts`): OS, both muons
p_T > 4.5 GeV & |η| < 2.4 & WP (NTuple stage), m_μμ ∈ (1.08, 2.9), pair p_T > 9,
**both muons `PassSingleMuFiducialGap(η, q)`, |η^pair| < 2.2**.

**ε_reco definition (NEW):**

    ε_reco(cell) = N[ truth pair in T  AND  both muons reco-matched  AND  pair passes the WP
                      AND the RECO pair passes the DATA signal selection incl. the gap cuts on RECO q·η and RECO η^pair ]
                 / N[ truth pair in T ]

    T = truth SIGNAL region, NO gap cut of any kind:
        single-b OS pair (`from_same_b`), truth p_T > 4.5 GeV & |truth η| < 2.4 both muons
        (NTuple-stage truth fiducial), truth m_μμ ∈ (1.08, 2.9), truth pair p_T > 9 GeV.

Binned in TRUTH (pair p_T, η^pair, ΔR) on the canonical coarse axes
(`ParamsSet::pair_pt_coarse_bins` × `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`
× `CommonEffcyConfig::dr_bins_edges_for_reco_effcy`), MC-event weighted, binomial errors —
exactly as before; only the denominator's cut list changes.

**Truth pairs outside the measured η^pair range.** The coarse η^pair axis spans [−2.2, 2.2] (it
tiles the data's surviving region, user decision 2026-09-07). A truth pair with |η^pair| ≥ 2.2
lands in the axis overflow of BOTH numerator and denominator and therefore enters no cell, no
ΔR-integrated fallback cell, and not the inclusive fallback (`TH3::Integral()` excludes
overflow; `Project3D` routes the η overflow to the 2D overflow — verified 2026-09-17). So:
- truth |η^pair| < 2.2, reco |η^pair| ≥ 2.2 → denominator yes, numerator no → **inefficiency**
  (the user's named edge case);
- truth |η^pair| ≥ 2.2, reco |η^pair| < 2.2 → outside the measured truth range; the data pair
  is still counted and corrected with its reco (edge) cell, i.e. it is **bin migration into
  the measured range**, the same class of effect as p_T or mass migration across the signal
  cuts, handled by the detector-response / unfolding step, not by ε_reco. Its size is measured
  from the MC and reported (Progress Log) so the user can judge whether it needs its own
  treatment.
- truth fails the signal cuts (mass / pair p_T) but reco passes → neither leg (unchanged
  convention).

### 3. Step-by-step method

(a) **RDF hist filling** (`RDFBasedHistFillingPythiaFullsim.cxx`, mirrored in the HIJING-overlay
per-centrality nodes): three columns with physics names —
`pass_signal_truth` = truth mass window && truth pair p_T (signal cuts ONLY);
`pass_gap_truth` = fiducial windows on truth q·η of both muons && truth |η^pair| < 2.2;
`pass_signal_reco` = both reco-matched && reco mass window && reco pair p_T && the gap cuts on
reco q·η of both muons && reco |η^pair| < 2.2 (= the data selection).
Filters: `_pass_signal_truth` (DENOMINATOR, truth signal cuts only),
`_pass_<wp>_and_signal_truth_and_reco` (NUMERATOR = WP && pass_signal_truth && pass_signal_reco),
`_pass_signal_truth_gapcut` (= pass_signal_truth && pass_gap_truth; the MC-vs-data SIGNAL
family only — see (d)). The RDF output carries a `TNamed` provenance marker
`pair_reco_eff_definition` naming this denominator definition.

(b) **Builder** `build_pp24_fullsim_pair_reco_eff.C`: ratio, fallbacks and census unchanged;
REFUSES an input file without the (a) marker (a stale file would silently deliver the old
fiducial definition — the textual-vs-inherited hazard of `mu_pt45_gap125_pairpt9_adoption.md`
R4); provenance string updated ("gap cuts on the RECO leg only; contains the gap acceptance").

(c) **Evaluator / crossx** (`PairRecoEffEvaluator.h`, `RDFBasedHistFillingPP.cxx`): mechanics
unchanged (evaluated at reco kinematics, fallbacks, floor); comments updated. **No ε_acc factor
anywhere.**

(d) **MC-vs-data SIGNAL family** (`plots/mc_data_compr/signal/`): an UNCORRECTED shape
comparison of data (reco gap cuts) with truth MC — it keeps the truth gap cuts, so it moves to
the `_pass_signal_truth_gapcut` filters (name already used with this meaning by
`RDFBasedHistFillingPowhegFullsim.cxx`). Content unchanged; plotter keys updated.

(e) **Rerun**: fullsim pp24 FULL RDF hist filling (stage 5) → stage-7 plots (unaffected in
content — they are the inclusive WP ratio — but refreshed) → builder → applied-map plots →
pp24 crossx pipeline (`pipeline_pp_crossx.sh SKIP_CONDOR=1`) → R_AA plots if the Pb+Pb inputs
are current.

### 4. Negative constraints

1. Do NOT apply any gap cut (single-muon window or |η^pair|) on the TRUTH leg of the
   denominator or the numerator of ε_reco.
2. Do NOT drop any gap cut from the RECO leg: the numerator's reco selection must stay
   bit-for-bit the data `signal_cuts`.
3. Do NOT apply a separate ε_acc anywhere (it would double count) — and do not quote one.
4. Do NOT touch the gap windows, the signal cuts, the NTuple truth fiducial, or any binning.
5. Do NOT change the MC trigger-efficiency selection (Steps 1–4 keep the gap cuts on both legs:
   ε_trig is conditional on reco-SELECTED pairs) nor the data tag-and-probe.
6. Do NOT extend the coarse η^pair axis to catch the out-of-range truth pairs (binning rule).

## Context

- Gap-cut inventory and history: `muon_gap_cuts_acceptance.md` (F7/F10/F11/F12/F17).
- Current selection values: `mu_pt45_gap125_pairpt9_adoption.md` §2.
- ε_reco product chain: `RDFBasedHistFillingPythiaFullsim.cxx` → `histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root` → `build_pp24_fullsim_pair_reco_eff.C` → `reco_eff/pair_reco_eff_pp24_full.root` → `PairRecoEffEvaluator` in `RDFBasedHistFillingPP.cxx` (`w_reco`).
- Stage-7 reco-eff plots (`plots/pp24_reco_effcy_plots/tight/{2d,ranged,signed,single_b_op_compr}`) are the INCLUSIVE WP ratio (`require_signal_cuts=false`), no signal or gap cut on either leg → unaffected by this change. Only the "applied" set (`tight/applied/`, from the ROOT product) changes.

## Scope

IN: pp24 Pythia fullsim FULL sample (the applied ε_reco); the HIJING-overlay per-centrality
mirror (code only — its test-sample outputs belong to a sibling session, rerun flagged, not
done here); builder, evaluator, crossx comments; MC-vs-data signal-family filter rename; docs
(`analysis_overview.md`, `signal_selection_change_impact.md`, `placeholder.md`,
`pythia_fullsim_pp.md`, `ParamsSet.h` STATUS block, `muon_gap_cuts_acceptance.md` closure note).
OUT: `RDFBasedHistFillingPowhegFullsim.cxx` legacy reco-eff path (obsolete, memory
`project_mc_sample_roles`; its `pass_signal_truth` still includes the truth gap — recorded, not
changed, so no stale-on-disk artefact is created); MC trig-eff; data T&P; Pb+Pb crossx (uses
the Run-2 placeholder reco-eff, untouched).

## Design Decisions

- **D1 — column/filter naming follows the user's terminology.** `pass_signal_truth` = signal
  cuts only; the gap cuts get their own column `pass_gap_truth`; the with-gap truth family is the
  suffixed `_pass_signal_truth_gapcut` (opt-in, suffixed variant — CLAUDE.md Binnings rule 4
  spirit). Consequence: the on-disk `*_pass_signal_truth*` histograms change meaning → a loud
  marker (`TNamed pair_reco_eff_definition`) is written by the filler and required by the builder.
- **D2 — out-of-range truth η^pair pairs are migration, not efficiency** (Physics Procedure §2).
  Forced by the binning rule (axis stays ±2.2) and consistent with how mass/p_T migration is
  already treated.
- **D3 — backup of the previous applied plots** (user instruction, mid-task): the applied set
  for the CURRENT selection had never been drawn (`tight/applied/` was absent; only the
  2026-08-24 4-GeV-era set existed under `old_4GeV_cut/`), so it is generated from the current
  `pair_reco_eff_pp24_full.root` (2026-09-10 14:10, fiducial definition) for both WPs and moved,
  together with a copy of that ROOT product, to ONE directory
  `plots/pp24_reco_effcy_plots/before_gap_acceptance/`.

## Implementation Plan

| # | Step | Physics § | Review | Status |
|---|---|---|---|---|
| 1 | Backup: draw the before-change applied plots (tight + medium) from the current ROOT product; move to `before_gap_acceptance/` with the ROOT file | D3 | — | DONE |
| 2 | `RDFBasedHistFillingPythiaFullsim.cxx`: columns + filters + provenance marker; overlay mirror | §3a | /review-analysis-code | DONE |
| 3 | Builder guard + provenance; evaluator + crossx comments; ParamsSet STATUS comment | §3b, §3c | /review-analysis-code | DONE |
| 4 | MC-vs-data signal family → `_pass_signal_truth_gapcut` (3 plotter files) | §3d | /review-analysis-code | DONE |
| 5 | Compile all touched classes (ACLiC); small test on a subset if feasible | — | — | DONE |
| 6 | Rerun stage 5 → 7, builder, applied plots; measure the migration-in fraction | §3e | /review-plot (covered by /review-analysis-code C1/C2 on the maps) | DONE |
| 7 | pp24 crossx pipeline (SKIP_CONDOR=1), R_AA if inputs current | §3e | /review-analysis-code C1/C2 on the crossx (it.2) | DONE (R_AA mixed-fiducial, R3) |
| 8 | Docs + INDEX; commit | — | /review-analysis-code APPROVED it.3 | DONE |

## Progress Log

- **Step 1 DONE (2026-09-17).** Before-change applied set drawn from the current
  `reco_eff/pair_reco_eff_pp24_full.root` (2026-09-10 14:10, fiducial definition) for Tight and
  Medium with `plot_pp24_fullsim_pair_reco_eff.cxx` (unchanged), then moved to
  `~/usatlasdata/pythia_fullsim_full_sample/plots/pp24_reco_effcy_plots/before_gap_acceptance/{tight,medium}/applied/`
  (5 PNGs + `pair_reco_eff_values.txt` each) together with a copy of the ROOT product and a
  `README.txt`. Before-change reference numbers (Tight): 143 of 288 3D cells populated; ε vs ΔR
  integrated 0.7755 / 0.7157 / 0.6340 / 0.3370 in the four ΔR bins.

- **Steps 2–4 code DONE (2026-09-17), compiling.** Files:
  `RDFBasedHistFilling/RDFBasedHistFillingPythia.h` (static `TruthGapCutExpr()` /
  `RecoGapCutExpr()` shared by pp24 + overlay), `RDFBasedHistFillingPythiaFullsim.cxx`
  (`pass_signal_truth` = signal cuts only; new `pass_gap_truth`; new
  `_pass_signal_truth_gapcut` filters for `_ss/_op/_single_b`; MC-vs-data signal family
  points at them; `TNamed pair_reco_eff_definition` written in `WriteOutputExtra`),
  `RDFBasedHistFillingPythiaFullsimOverlay.cxx` (per-centrality mirror, no `_gapcut` family),
  NEW `Utilities/PairRecoEffDefinition.h` (Key/Value of the marker),
  `plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C` (refuses input without the
  marker → return 4; copies it into the product; provenance text),
  `Utilities/PairRecoEffEvaluator.h` (refuses product without the marker; definition comment),
  `RDFBasedHistFillingPP.cxx` (comments), `ParamsSet.h` (STATUS comments),
  `plotting_codes/mc_data_compr/{PlotMCDataComprBaseClass.h,plot_mc_data_pair_pt_in_eta.cxx,
  plot_mc_data_compr_signal.cxx,McDataComprConfig.h}` (`_gapcut` keys),
  `plot_pp24_fullsim_pair_reco_eff.cxx` (values.txt header). Docs: `analysis_overview.md`
  §2/§4b, `signal_selection_change_impact.md`, `placeholder.md`, `pythia_fullsim_pp.md`,
  `muon_gap_cuts_acceptance.md` F19, `INDEX.md`.

- **Steps 5–7 DONE (2026-09-17).** Compile OK (3 classes + builder + plotter). Builder refused
  the pre-change file (return 4). `pipeline_pythia_fullsim_pp.sh full` with `SKIP_NTP=1`:
  stages 4–9 DONE 20:43 (`pipelines/logs/fullsim_pp_gap_acc_20260917.log`). Product rebuilt
  (20:44), applied plots drawn for Tight + Medium. `pipeline_pp_crossx.sh` (`SKIP_CONDOR=1`)
  completed 20:48 (`pipelines/logs/pp_crossx_gap_acc_20260917.log`; evaluator: 288 cells, 145
  empty → dR-integrated fallback, 0 inclusive fallback, 0 floored; MC-vs-data stage 8 read the
  `_gapcut` keys). `RAA_plotting.cxx` rerun 21:03 → `plots/single_b_analysis/RAA/raa_{pair_pt,
  pair_eta,ctr}_pbpb23_24_25_26_combined_pp24_2mu4.png` (the macro's default now combines the
  four Pb+Pb years present on disk). **See R3 — R_AA is now MIXED-FIDUCIAL.**
- **Step 8 review DONE:** `/review-analysis-code` APPROVED at iteration 3
  (`.claude/logs/review-analysis-code-20260917-204616-reco-eff-gap-acceptance.md`): all numbers
  independently reproduced; iteration-1 WARNING (truth region must name |η^pair| < 2.2 as the
  measured range) and iteration-2 WARNING (same statement in `signal_selection_change_impact.md`)
  fixed; C1/C2 pass on the reco-eff maps and on the regenerated pp24 crossx; C3 UNVERIFIED (no
  comparable Run 2 pair-with-gap number in the KB).

## Results & Observations

### R1 — New ε_reco vs the fiducial one (2026-09-17)

Source: refilled `histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root` (stage 5,
2026-09-17 20:39; previous file backed up by the pipeline to
`backup/histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.bak_20260917_203826.root`),
carrying the `pair_reco_eff_definition` marker. Product rebuilt →
`reco_eff/pair_reco_eff_pp24_full.root` (the builder had REFUSED the pre-change file: return 4,
old product untouched — guard verified).

| | Tight | Medium |
|---|---|---|
| inclusive ε_reco, OLD (gap cuts on both legs) | 0.7157 | 0.7780 |
| inclusive ε_reco, NEW (gap cuts on reco leg only) | **0.5850** | **0.6360** |
| ratio NEW/OLD = pair-level gap acceptance folded in | 0.817 | 0.817 |
| 3D cells with no measure | 145 / 288 | 145 / 288 |
| 2D (pT,η) cells empty | 0 / 72 | 0 / 72 |

The ratio is WP-independent to 3 decimals, as it must be (the gap loss is a truth-level
acceptance effect, the WP acts only on reconstructed muons).

**Migration across |η^pair| = 2.2 (Physics Procedure §2, D2), measured** (scratch macro
`migration_in.C`, MC-weighted): truth |η^pair| ≥ 2.2 pairs in the η overflow of the numerator =
**0.02 % of the in-range numerator (0.24 % of the two edge panels)** for both WPs. Negligible —
no dedicated treatment needed; it is part of the detector-response step in any case. The truth
DENOMINATOR weight in the η overflow is **3.35 % relative to the in-range denominator (0.1919 /
5.723; 3.2 % of the whole truth signal region)** — those pairs are outside the measured range,
so |η^pair| < 2.2 is a truth fiducial edge of the measurement (single-muon windows
acceptance-corrected, pair-level edge not; stated in `analysis_overview.md` §2).

### R2 — Unchanged by construction
Stage-7 reco-eff plots (`tight/{2d,ranged,signed,single_b_op_compr}`) are the inclusive WP
ratio and do not read the signal legs; detector-response, single-muon reco-eff and the kn
tables are untouched. The MC-vs-data SIGNAL family (`_pass_signal_truth_gapcut`) has the same
content as the old `_pass_signal_truth` family under a new name.

### R3 — ★ R_AA is MIXED-FIDUCIAL until the Pb+Pb pair reco-eff adopts the same definition

pp24 is now corrected for the gap acceptance (per coarse η^pair panel: ×1/0.45 in |η^pair| ∈
[1,1.5], ×1/0.48 in [−0.5,0.5], ×1/0.79 elsewhere). Pb+Pb still uses the Run-2 single-muon
PLACEHOLDER reco-eff (`reco_eff_placeholder_run2.md`), which contains NO gap acceptance, while
the Pb+Pb signal selection applies the same gap cuts (D1) — so the Pb+Pb crossx is still a
gap-FIDUCIAL one. Consequence, visible in `raa_pair_eta_pbpb23_24_25_26_combined_pp24_2mu4.png`:
R_AA vs η^pair now carries DIPS of ≈0.6 (|η^pair| ∈ [1,1.5]) and ≈0.6–0.7 ([−0.5,0.5]) relative
to the neighbouring panels — the pp acceptance correction with no Pb+Pb counterpart — and the
integrated R_AA is biased LOW by ≈0.82. This is NOT a bug of this change; it is the pp/Pb+Pb
asymmetry of the reco-eff inputs, and it resolves only when the Pb+Pb pair reco-eff (HIJING
overlay, per centrality; code mirror already in place) is built with the same definition. Until
then R_AA must not be quoted (it was already PROVISIONAL per `mu_pt45_gap125_pairpt9_adoption.md`).

## Remaining Work

1. **Pb+Pb pair reco-eff with the same definition** — the overlay per-centrality mirror is coded
   and compiles, but its products are test-sample only and the Pb+Pb crossx reads the Run-2
   placeholder. Needs the FULL overlay production (user/production item) and then the Pb+Pb
   analog of the builder/evaluator. Until then R_AA is mixed-fiducial (R3).
2. **Overlay test-sample histograms are STALE w.r.t. the code** (pbpb23 test 2026-09-09, pbpb24
   test 2026-09-16, both filled with the old `pass_signal_truth` = signal+gap definition and no
   marker). Not refilled here — the pbpb24 test dir belongs to the active
   `hijing_overlay_pbpb24_test_sample_skim.md` session. Refill when convenient; nothing
   downstream consumes them today.
3. `RDFBasedHistFillingPowhegFullsim.cxx` legacy reco-eff path keeps `pass_signal_truth` =
   signal + truth gap (documented divergence; obsolete path). Align if that path is ever revived.
4. Note sections (`IntNotes/tex/reconstruction_efficiency.tex`) describe the fiducial definition
   — owned by the active `intnote_trig_reco_eff_update_2026_09.md`; flagged there is NOT done
   here (cross-session), so that doc must pick up R1/R3 + `analysis_overview.md` §4b.

## Latest Stage

**DONE 2026-09-17 (pp24 side complete; Pb+Pb side pending — Remaining Work 1–2).** Definition
change implemented, marker-guarded, reviewed (APPROVED it.3), pp24 fullsim RDF + product +
applied plots + pp24 crossx + MC-vs-data + R_AA regenerated. Before-change applied plots and
product kept in `plots/pp24_reco_effcy_plots/before_gap_acceptance/`. Headline: inclusive
ε_reco 0.7157 → 0.5850 (Tight), 0.7780 → 0.6360 (Medium), ratio 0.817 both WPs. **R_AA is now
mixed-fiducial (R3) until the Pb+Pb pair reco-eff adopts the same definition.**
