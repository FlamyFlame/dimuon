# MC-vs-data comparison plots: split into SIGNAL-region and GENERIC families, with matched cuts

Mode: **IMPLEMENTATION**. Opened 2026-08-25.

## Objective

Rebuild the pp24 MC-vs-data comparison plot set
(`/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/`) so that

1. it is split into two clearly separated families in **separate subdirectories** —
   `signal/` (data-like **signal-region** cuts, the urgent deliverable) and `generic/`
   (no signal-region cuts);
2. **every sample in a family carries the SAME cuts** — Pythia (and POWHEG in `generic/`)
   must mirror the data selection, in truth quantities where the sample is truth, plus
   `from_same_b` for the single-b signal definition;
3. Pythia comes from a **centrally produced FULL sample**, never the private sample;
4. the known defects of the current set are fixed (unity variants removed, data
   "jacobian corrected" histograms, POWHEG multiple counting, stale turn-on file), and
   orphaned legacy plot directories are deleted.

Parent context: `pp24_crossx_rerun_2026_08.md` (DONE 2026-08-18) produced the current set and
left OPEN items **O1** (private-Pythia normalization is a tuned empirical factor) and **O2**
(eps_reco extrapolated over the generic domain). This doc executes O1 and documents O2.

## Autonomy Contract (DONE 2026-08-25 — all six Done items met)
- Mandate: run autonomously to DONE; do NOT pause to confirm progress. Finishing a
  plan, a passing small test, or one pipeline stage is NOT a stopping point.
- Done =
  1. This doc carries a complete, numbered execution plan covering ALL ten user items,
     with physics justification where the change is physics.
  2. **SIGNAL family (PRIORITY)** produced at `plots/mc_data_compr/signal/`: pair_pt,
     pair_pt_in_eta_subplots, DR_zoomin, Deta_zoomin, Dphi_zoomin, minv_zoomin, pair_eta
     (-2.4..2.4) — data and Pythia on IDENTICAL cuts (data: OS + Tight + 1.08<m<2.9 +
     pair pT>8 + fiducial gap cut on both muons, efficiency-corrected; Pythia: central FULL
     sample, truth kinematics, `from_same_b` + the same mass / pair-pT / gap cuts,
     AMI-normalized). No unity variants, no jacobian variants.
  3. **GENERIC family** produced at `plots/mc_data_compr/generic/` with Pythia (and POWHEG)
     carrying the same non-signal cuts as the data, the data jacobian histograms actually
     jacobian-corrected, the POWHEG multiple-counting removed, and a `README.md` in that
     directory warning that the corrections are extrapolated outside their measurement region.
  4. All histogram inputs regenerated in ONE grouped RDF rerun per sample, using the CURRENT
     single-muon turn-on fit file.
  5. Unity-variant plots removed from the code and from disk; orphaned subdirectories of
     `plots/mc_data_compr/` deleted.
  6. Docs updated (this doc, INDEX.md, `signal_selection_change_impact.md` if the map changed)
     and the work committed.
- Stop-and-ask = ANY physics-results-bending ambiguity (no fixed list; use judgment;
  when unsure whether an ambiguity is blocking, treat it as blocking → AskUserQuestion).

## Physics Procedure

### 1. Motivation

A data/MC comparison is only interpretable if the two curves describe the **same phase-space
region** and the **same object**. The current set violates this in four ways at once:

* the data are efficiency-corrected inside a fiducial region defined by the detector-gap cut,
  while the Pythia truth curve has no gap cut at all (an `eps_acc`-sized, ~16 % offset —
  `pp24_crossx_rerun_2026_08.md` M5);
* the data are corrected with efficiencies measured ONLY inside the single-b signal region,
  but the generic histograms have no signal-region cut, so the correction is extrapolated
  (O2);
* the Pythia curve comes from a PRIVATE sample whose weight never touches AMI, normalized by
  a hand-tuned 1e6 factor (O1), so the absolute level carries no physics;
* the POWHEG curve sums three OVERLAPPING classification axes, over-counting by ~2.88x.

Splitting into a signal family and a generic family is what makes the first two tractable:
inside the signal region every correction is applied where it was measured, so the comparison
is a genuine (fiducial, reco-level, background-unsubtracted) cross-section comparison. The
generic family keeps the shape information over the full phase space but is explicitly labelled
as extrapolated.

### 2. Top-level statement

**SIGNAL family.** Per observable X in {pair pT, pair eta, dR, dphi, deta, m_uu}:

```
DATA:    dsigma/dX = (1/L_int) * SUM_{OS pairs in signal region} 1/[eps_trig^pair * eps_reco^pair]
PYTHIA:  dsigma/dX = SUM_{truth single-b OS pairs in the SAME region, truth kinematics} w_AMI * 1000
```

with `L_int = 400.412 pb^-1` (pp24, 2mu4), the efficiency factors exactly as in
`pp24_crossx_rerun_2026_08.md` §2/§3, and the nb->pb factor 1000 because the AMI
`crossSection` is in nb. NO unfolding, NO eps_acc, NO background subtraction — the data curve
is a fiducial reco-level OS yield and the Pythia curve is pure truth single-b; both facts must
be quoted with any number read off these plots.

**GENERIC family.** Same weights, but with the signal-region cuts (mass window, pair-pT
threshold, `from_same_b`) removed on BOTH sides. The fiducial gap cut stays on both sides — on
the data side it is mandatory (a muon at q*eta in [2.30,2.40) has no fitted turn-on and the
evaluator throws), and on the MC side it is required for the two curves to describe the same
region. Absolute normalization here is NOT a cross-section comparison for the data side,
because eps_reco and eps_dR are evaluated outside the domain in which they were measured
(see §4 and the directory README).

### 3. Step-by-step method

**(a) The data signal region** — the single source of truth is the pp24 crossx selection,
read from code, not retyped: opposite sign, Tight working point, `1.08 < m_uu < 2.9` GeV,
`pair pT > 8` GeV, and BOTH muons outside every window of
`ParamsSet::single_mu_fiducial_gap_cuts`. No dR cut (removed 2026-06-22).
[TO BE CONFIRMED FROM CODE — see Progress Log Step 1.]

**(b) The Pythia mirror.** The same cuts in TRUTH quantities, plus `from_same_b` (the single-b
signal definition; the data has no truth information and its OS yield therefore contains
gluon-splitting and combinatorial background — a deliberate, documented mismatch that CANNOT be
removed without the template fit). The gap cut is applied to TRUTH `q*eta`, mirroring the
`eps_reco` denominator (`pp24_crossx_rerun_2026_08.md` §3d), so that the chain
(truth MC region) == (data fiducial region after correction) closes.

**(c) Binning.** All 1D observables use the binning already registered for the data side in
`var1D_pp.json`; pair pT uses `ParamsSet::pT_bins_150`; pair eta uses the 44 uniform bins of
the crossx 2D histogram, and the 9-panel view uses
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`. **No binning is invented and none is
retyped** — the MC histograms are booked from the SAME binning entries as the data ones so that
a bin-edge divergence is impossible by construction.

**(d) Normalization.** Data: `1/L_int` via `PPBaseClass::GetCrossxFactor(24,"2mu4")`.
Pythia: the AMI weight already carried by the central full sample, times 1000 (nb->pb).
The private sample's tuned `1e6` factor is DELETED with the private sample.

### 4. Negative constraints (things the code must NOT do)

1. Do **NOT** apply a different cut set to MC than to data within one family. The whole point
   of the split is that each family is one region.
2. Do **NOT** use the private Pythia sample anywhere.
3. Do **NOT** re-invent or retype any binning (pair pT, pair eta, q*eta, gap windows, the 1D
   axes) — read them from `ParamsSet` / `CommonEffcyConfig` / the json binning configs.
4. Do **NOT** apply `eps_acc`; the result stays fiducial, exactly as the cross-section is.
5. Do **NOT** present the generic family's absolute level as a cross-section: eps_reco is
   measured only inside the signal region and is CLAMPED at the axis edges outside it, and
   eps_dR is exactly 1 outside its (pair pT, pair eta) grid and above dR = 1. The dR>1 part is
   physically defensible — the dR correction originates in resolution/RoI effects that only
   act on nearby pairs — but the pair-pT and mass extrapolation is not. This warning goes in
   `generic/README.md`.
6. Do **NOT** sum overlapping POWHEG classification axes (mechanism / origin / flavour are
   three views of the SAME sample, not a partition).

## Context — observations inherited from the 2026-08-25 read-only investigation

The set as it stood on 2026-08-18 was traced end to end. Findings, all verified against code
or against the ROOT files themselves:

**F1. The data `*_jacobian_corrected` histograms are NOT jacobian corrected.**
`RDFBasedHistFillingData.cxx:211` (and `:238` for the `_wgapcut` twin) calls the
`(suffix, df, weight_col, ...)` overload with `generic_weight_col = w_reco_trig`, bypassing the
`weight_specifier_to_column_map` lookup; no 1/dR weight column exists on the data side at all.
Verified in ROOT: `h_DR_op` and `h_DR_op_jacobian_corrected` are bit-identical (integral
1.2727e7; bin 5 = 264460 in both). The Pythia one IS 1/dR weighted (`weight_over_dr`).
=> both `*_jacobian_corrected` PNGs overlaid a 1/dR-weighted Pythia on an UNweighted data and
POWHEG curve.

**F2. The POWHEG curve is multiply counted (~2.88x) and shape-distorted.**
`PlotMCDataComprBaseClass.h::powheg_mechanisms` sums 10 histograms that are THREE overlapping
axes, not a partition: mechanism {gg, qg, single_g, qq, others}, origin {from_same_b, single_b
— numerically IDENTICAL, 4346.05 each}, flavour {bb, cc, other_flavors}. Verified on
`h_DR_sign2_*`: true partition total 9864.3 vs sum-of-10 = 28420.7 -> 2.88x, with single-b
entering twice extra. That POWHEG file is dated Mar 2025 and predates the current producer
naming (`_origin_binned_*` / `_flavor_binned_*` in `RDFBasedHistFillingPowhegTruth.cxx`).

**F3. Working-point consistency in the generic family.** The generic histograms are booked
before the crossx `pair_pass_tight` filter. Whether that is a real WP mismatch depends on
whether the TIGHT WP is already required at the NTUPLE stage — to be settled from code in
Step 1 (the user's recollection is that the pp crossx pipeline moved the Tight requirement to
the ntuple level, in which case there is no mismatch and no RDF-stage filter is needed).

**F4. Corrections extrapolated outside their measurement region (generic family only).**
eps_reco was measured only inside the signal region and the evaluator CLAMPS pair pT / pair eta
/ dR to the edge cells; eps_dR returns 1 outside its grid and above dR = 1. Census from the
2026-08-18 production run: 73.6 % of evaluations clamped in dR, 52.5 % in pair pT, 21.8 % in
pair eta, 12.2 % fell back to the dR-integrated map. This is `pp24_crossx_rerun_2026_08.md` O2.

**F5. The turn-on fit file moved after the plots were made.**
`pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit.root` is dated
2026-08-24 17:14, six days AFTER the 2026-08-18 hist filling, so the plots on disk use the
previous fit. Any rerun picks up the new one automatically.

**F6. The 1D family's Pythia normalization is meaningless in absolute terms** — private
sample, `eventWeight/njobs` only, no AMI, no isospin, no sigma, times a hand `1e6` whose second
factor of 1e3 is an acknowledged truth-combine under-normalization bug. This is
`pp24_crossx_rerun_2026_08.md` O1, and the user's instruction to move to the central full
sample resolves it.

**F7. The two families read DIFFERENT samples and different histogram types today.**
Family A (10 shape overlays, `plot_mc_data_compr.cxx` + `PlotMCDataComprBaseClass.{h,c}`) reads
the private Pythia truth sample and the generic 1D data histograms `h_<kin>_{ss,op}`.
Family B (2 cross-section plots, `plot_mc_data_pair_pt_in_eta.cxx`) reads the central FULL
Pythia truth sample (`h2d_sig_accept_num_pt_150_eta`) and the crossx 2D data histogram
`h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts`. Input file dates: data 2026-08-18,
central Pythia 2026-06-23, private Pythia 2026-03-06, POWHEG 2025-03-03.

## Implementation Plan

Ordering principle (user): (i) SIGNAL family first and physically correct; (ii) group the
time-costly work so each sample's RDF stage is rerun ONCE.

**Step 1 — ground truth from code.** DONE (Progress Log 1a-1d).

**Step 2 — shared binning registration** (orchestrator-owned; both other steps depend on it).
Add `hist_binning_map["pair_eta_crossx"]` = 44 uniform bins on [-2.4,2.4], built once in
`RDFBasedHistFillingBaseClass.cxx::BuildHistBinningMap`, and replace the four retyped literals
with it. Per Physics Procedure §3c / D2.

**Step 3 — data-side RDF** (`RDFBasedHistFillingPP.cxx`, `RDFBasedHistFillingData.cxx/.h`,
`var1D_pp.json`). Per §3a and D3:
  a. hoist the `pair_pass_tight` filter above the generic booking (D3);
  b. add the SIGNAL-REGION 1D histogram set — OS and SS, filter = `signal_cuts`, weight =
     `crossx_weight_trig_corr` (so the histogram IS d(sigma) in pb) — for
     `DR_zoomin`, `Dphi_zoomin`, `Deta_zoomin`, `minv_zoomin`, `pair_eta_crossx`,
     `pair_pt_log_150`, every binning READ from `var1D_pp.json` / `hist_binning_map`;
  c. fix F1: define a real `w_reco_trig_over_dr` column and route the `_jacobian_corrected`
     booking through it.

**Step 4 — MC-side RDF** (`RDFBasedHistFillingPythiaFullsim.cxx`, `var1D_pythia_fullsim.json`).
Per §3b and D1/D2: add the SAME six observables, in TRUTH quantities, on the DATA binning
(`*_ppbin` / `truth_pair_eta_crossx` / `truth_pair_pt_log_150`), under
  * `_single_b_pass_signal_truth` (SIGNAL family — the cut set is already exactly the data's), and
  * a new truth-gap-cut filter on `_op` and `_ss` (GENERIC family: gap cut only, no mass window,
    no pair-pT threshold, no `from_same_b`, mirroring the data generic selection).

**Step 5 — generic-family code fixes** (`plotting_codes/mc_data_compr/`): remove every `_unity`
variant, fix the POWHEG hist list to the flavour axis alone (D6), drop the jacobian panels from
the signal family, delete the `TRUTH_COMBINE_RENORM` hand factor (D1).

**Step 6 — ONE grouped rerun per sample**: pp24 data RDF crossx hist filling (picks up the
CURRENT `single_mu_effcy_pT_fit.root`, F5, and every Step-3 change at once) and the pp24 fullsim
RDF stage (~42 s). No NTuple rerun is needed on either side — every column already exists.

**Step 7 — plotting**: output-directory mechanism; `signal/` (pair_pt, pair_pt_in_eta_subplots,
DR_zoomin, Deta_zoomin, Dphi_zoomin, minv_zoomin, pair_eta) and `generic/` (the surviving shape
overlays + POWHEG). Pipeline Stage 8 extended to run both macros.

**Step 8 — housekeeping**: delete the 8 orphaned subdirectories; write `generic/README.md` (§4
constraint 5 + D4); update `analysis_overview.md` / `signal_selection_change_impact.md` if the
map moved; INDEX.md; commit.

## Design Decisions

**D1. The Pythia partner is the pp24-condition FULLSIM FULL sample, for BOTH families.**
File: `/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root`
(producer `RDFBasedHistFillingPythiaFullsim.cxx` <- `pipelines/pipeline_pythia_fullsim_pp.sh full`).
Rejected alternative: the Pythia TRUTH full sample (`pythia_truth_full_sample/pythia_5p36TeV/`).

*Physics reason (decisive).* The truth full sample's weight sets `nominal_beam_ratio` to the Pb
isospin mixture **4:6:6:9 unconditionally** (`PythiaAlgCoreT.c:344-347`) and the production runs
all four beams, so its absolute sigma is a **Pb isospin-averaged NN cross-section, not a pp
cross-section** — the wrong object to put beside pp24 data. The fullsim pp24 sample is pp beam
only with isospin weight 1 (`PythiaAlgCoreT.c:581-590`). Measured size of the difference in the
signal region: truth-sample 9.408 vs fullsim 7.844 (same units) = **16.6 %**, i.e. the same order
as the whole data/MC discrepancy being studied. The existing `_pp_only` truth NTuples do NOT fix
it (they drop the pn/np/nn chains but the surviving pp beam keeps its 4/25 share, so they are
6.25x low) and have no RDF histogram output at all.

*Three supporting reasons.* (i) The fullsim class ALREADY implements the data signal region
bit-for-bit on truth: `pass_signal_truth = truth_minv > 1.08 && truth_minv < 2.9 &&
truth_pair_pt > 8 && FiducialGapCut(truth q*eta, both muons)`
(`RDFBasedHistFillingPythiaFullsim.cxx:178-187`), identical to the data `signal_cuts`
(`RDFBasedHistFillingPP.cxx:531-534`) — nothing has to be re-derived, and the site that would
have to be changed on the truth side (`RDFBasedHistFillingPythiaTruth.cxx:326-328`) is the one
that produces the low-mass template-fit `_sigsel` templates, i.e. it carries a blast radius into
`low_mass_dimuon_template_fit.md`. (ii) It carries truth AND reco kinematics plus `from_same_b`
(`MuonPairPythiaFullSimWTruth`, `MuonPairPythia.h:76-84`), and its pairs are TRUTH-SEEDED
(`PythiaFullSimExtras.c:371,419-425`), so truth distributions are not reco-biased. (iii) It is
the most recently maintained sample (2026-08-18) and its RDF stage reruns in ~42 s.

*Accepted cost, to be quoted with the plots.* The fullsim has ~3.8x / 6.8x fewer events than the
truth sample in the top two pT-hat slices, so the 50-150 GeV pair-pT tail is statistically
weaker. It shows as error bars, not as a bias.

*Bonus.* The fullsim `weight` is a genuine AMI weight in nb, so the `TRUTH_COMBINE_RENORM = 1e3`
hand factor of `plot_mc_data_compr.cxx:130` is **deleted**; only the physical nb->pb x1000
survives. This closes `pp24_crossx_rerun_2026_08.md` **O1**.

**D2. Binning: the comparison uses the DATA's registered binning on both sides; nothing existing
is changed.** Per `.claude/CLAUDE.md` Binnings rule 1 the values are READ from their registry and
never retyped, and per rule 4 the MC-side comparison histograms are **opt-in, suffixed variants**
(`*_ppbin`) that sit beside — and do not replace — the MC's native axes, which keep serving the
MC-only reco-efficiency and detector-response plots. Concretely the MC and data 1D axes disagree
today (fullsim `truth_dr_zoomin` 20x[0,1] vs data `DR_zoomin` 20x[0,0.8]; `dphi/deta_zoomin`
+-1.0 vs +-0.8; `pair_eta` 48 vs 44 bins), which is why a ratio pad was impossible and the
overlays were only density comparisons.
  * pair pT -> `ParamsSet::pT_bins_150` (already shared).
  * pair eta -> the crossx fine axis, 44 uniform bins on [-2.4, 2.4]. This value is currently a
    RETYPED literal in four places (`RDFBasedHistFillingPP.cxx:693,757,760`,
    `RDFBasedHistFillingPbPb.cxx:1139,1294`). It is promoted to ONE named binning,
    `hist_binning_map["pair_eta_crossx"]`, built once in the base class and referenced by name
    from the crossx booking and from both json var lists. Same numbers, one definition.
  * dR / dphi / deta / minv zoom-ins -> the `var1D_pp.json` entries (`DR_zoomin` 20x(0,0.8),
    `Dphi_zoomin` 20x(-0.8,0.8), `Deta_zoomin` 20x(-0.8,0.8), `minv_zoomin` 40x(0,3.0)).

**D3. The Tight WP filter is hoisted above the generic booking.** It currently sits inside
`FillHistogramsCrossx` (`PP.cxx:543-546`) while the generic histograms are booked first
(`Data.cxx:151-158`), so the generic family is a MEDIUM-WP yield divided by TIGHT efficiencies.
Moving the filter to where the fiducial gap cut is applied (`PP.cxx:470-478`) makes both families
Tight; the crossx call then becomes an idempotent no-op. This is a correction, not a convention
change: `eps_reco` (0.6620) and the turn-on fits in force are the Tight ones.

**D4. The generic family's resonance-veto mismatch is documented, not fixed.** The data OS tree is
resonance-vetoed at the ntuple stage (a muon in an OS pair with `minv > 60` or minv in
`ParamsSet::minv_cuts` tags BOTH its muons, and every pair containing a tagged muon is dropped);
the MC file is `_no_data_resonance_cuts` and a `_with_data_resonance_cuts` variant would need a
full MC NTuple Condor rerun of a 4.5 GB tree. In the SIGNAL region the veto is near-inert (it can
only remove a signal pair through a THIRD muon that formed a vetoed pair with one of its legs);
over the generic mass range it is a real mismatch. Recorded in `generic/README.md`.

**D5. Directory split.** `plots/mc_data_compr/signal/` and `plots/mc_data_compr/generic/`.
`PlotMCDataComprBaseClass` gains the output-directory mechanism it never had (today every path is
a string literal inside the `SaveAs` calls).

**D6. POWHEG appears in `generic/` only, on ONE classification axis.** The current
`powheg_mechanisms` list sums mechanism {gg,qg,single_g,qq,others} + origin {from_same_b,
single_b} + flavour {bb,cc,other_flavors} — three overlapping views of one sample, hence the
2.88x. The plot's own legend distinguishes bb from cc, so the **flavour** axis is the one kept.

**D7. Colour convention of the SIGNAL family (user instruction, 2026-08-25).**
`pp data` = **red**, `Pythia` = **black** for the curve that is the data's counterpart on each
pad (single-b on the OS pad, all-SS on the SS pad); the secondary `Pythia, all OS` reference
stays **blue**, deliberately not black, so it can never be read as the single-b signal.
Implemented as a signal-family-local palette (`kSignalDataColor` / `kSignalMcColor` /
`kSignalMcRefColor` in `PlotMCDataComprBaseClass.h`) rather than by editing the shared
`colors[]` array, so the GENERIC family's palette is unchanged.

**D8. POWHEG is brought onto the data's cuts ADDITIVELY (Step 9).** `RDFBasedHistFillingPowhegTruth.cxx:158`
applies the one-sided `q*eta < 2.2` and is the producer of the NLO template-fit input, so it is
NOT changed. A new gap-cut generic filter is added beside it, plus the six comparison variables
on the shared named axes — the same pattern used for the fullsim. This closes both POWHEG
mismatches (the missing gap cut and the 100-vs-40-bin dR axis) without touching
`low_mass_dimuon_template_fit.md`.


## Progress Log

- 2026-08-25 Step 0: doc created; Autonomy Contract pinned; the 2026-08-25 read-only
  investigation's findings F1-F7 recorded above as the starting context (they were previously
  only in conversation).


- 2026-08-25 Step 1a (scout: orphaned plot directories). All 12 top-level PNGs ARE producible
  (10 by `plot_mc_data_compr.cxx()`, 2 by `plot_mc_data_pair_pt_in_eta.cxx+()`) — no orphan at
  top level. **Eight subdirectories, 117 PNGs, are ORPHANED** (no code path writes there today):

  | dir | files | newest | why orphaned |
  |---|---|---|---|
  | `backup/` | 12 | 2023-03-07 | hand-made copy |
  | `backup_no_effcy_corr_w_mu4_mu4noL1/` | 32 | 2026-06-10 | hand-made bulk `cp` (all 32 within one second); holds the removed `_ratio_to_pp17` modes |
  | `dphi/` | 7 | 2023-03-02 | legacy `dphi_plots/*.c` write top-level relative paths; their input dir `usatlasdata/athena/runMCV2/` no longer exists |
  | `hist_2D/` | 25 | 2025-01-30 | see below |
  | `hist_2D_interval_projections/` | 6 | 2025-03-20 | see below |
  | `hist_2D_interval_projections/plots_with_4_to_8_GeV/` | 5 | 2025-03-20 | string appears nowhere in the repo |
  | `noDPHI_div/` | 12 | 2023-03-07 | producer `others/plot_powheg_data_compr_noDPHI_div.c` deleted in `6b952b7` |
  | `powheg_data_compr/` | 17 | 2023-03-23 | producer `others/plot_powheg_data_compr.c` deleted in `6b952b7` |
  | `with_gapcut/` | 7 | 2023-03-22 | producer `others/plot_mc_data_compr_w_gapcut.c` deleted in `6b952b7`; `with_gapcut` has 0 grep hits |

  `plot_mc_data_2D_hists_and_1D_proj.cxx` still CONTAINS the `hist_2D/` and
  `hist_2D_interval_projections/` path literals but **cannot run today**, for three independent
  reasons: (1) `hist_2D/` is behind `draw_2d_hist`, which `initialize()` forces false whenever
  `proj_intervals` is non-empty, and both live callers pass non-empty; (2) the histogram name it
  builds, `h_minv_pair_pt_zoomin_near_sign1`, exists in neither input file (the pp24 file has 10
  TH2, none matching; Pythia uses `h_minv_10GeV_vs_pair_pt_sign1_<mech>`) -> hard throw at line 87;
  (3) its MC entry sets `pythia_with_data_resonance_cuts = true`, asking for
  `histograms_pythia_combined_with_data_resonance_cuts.root`, which does not exist. So the two
  directories are orphaned in practice even though a macro names them.

  Also established: `PlotMCDataComprBaseClass` has **no output-directory mechanism at all** — every
  path is a string literal inside the `SaveAs` calls; `pythia_with_data_resonance_cuts_dir`
  (`PlotMCDataComprBaseClass.c:6`) is a dead member referenced nowhere. Adding the `signal/` and
  `generic/` subdirectories therefore requires introducing that mechanism (Step 7).
  Pipeline coupling: `pipeline_pp_crossx.sh` Stage 8 runs ONLY `plot_mc_data_compr.cxx()`;
  `plot_mc_data_pair_pt_in_eta.cxx` is manual-only and must be added.

- 2026-08-25 Step 1b (scout: the pp24 data selection, from code). **Answers the user's question
  "remind me what signal cuts we are applying now".**

  **NTuple stage** (`NTupleProcessingCode/DimuonDataAlgCoreT.c`, loop `:628-827`; `PPExtras`
  defines no hooks so every pp cut is in the core), in code order:
  1. `:672` -> `:501-512` HLT_2mu4 fired AND the pair itself matched, `dimuon_b_2mu4_mindR`
     (dR_trig < `mindR_trig` = 0.02). Governed by `trigger_mode` (`.h:335`, = 3 from
     `run_pp_24_nominal.sh`).
  2. `:581` combined muon (`quality & 1`), both.
  3. `:583-587` **`requireTight ? (quality&16 : Tight) : (quality&8 : Medium)`** — and
     **`requireTight` defaults to FALSE** (`.h:340`); `run_pp_24_nominal.sh` sets only
     `trigger_mode = 3`. Corroborated by the file name: `:396`
     `tight_suffix = requireTight ? "_tight" : ""` and the input on disk is
     `muon_pairs_pp_2024_2mu4_mindR_0_02.root`, no `_tight`.
  4. `:589` IDCuts (`&32`), `:590` MuonCuts (`&256`), both muons.
  5. `:593` `|eta| < 2.4`, `:596` `pT > 4` GeV, both muons.
  6. `:601` **one-sided** `dP/P < 0.12` (`pms.deltaP_overP_thrsh`), no `fabs`.
  7. `:605-608` `|d0| < 2 mm`, `|z0 sin(theta)| < 2 mm`, both muons.
  8. `:612-614` track-charge consistency — `turn_on_track_charge` is OFF.
  9. `:759-762` photo-production veto is `if (isPbPb)` ONLY — not applied to pp.
  10. `:782-796` resonance veto v1 (`resonance_cut_mode = 1`): a muon is tagged if it belongs to
      an OS pair with `minv > 60` or minv inside `pms.minv_cuts` =
      {0,1.06},{2.9,3.3},{3.55,3.8},{9.08,10.5}; any pair containing a tagged muon is dropped.
  11. `:820` sign split: **sign1 = SS, sign2 = OS**.
  `pair_pass_tight = (m1.quality & m2.quality & 16)` is written unconditionally as a BRANCH
  (`:677`). **No mass window, no pair-pT threshold, no dR cut and no gap cut at the ntuple stage.**

  **RDF crossx stage** (`RDFBasedHistFillingPP.cxx`), the SIGNAL REGION:
  1. `:543-546` `if (isTight) Filter("pair_pass_tight")` — `isTight` defaults **true**
     (`RDFBasedHistFillingData.h:186`).
  2. the fiducial gap cut, inherited from the generic pass which ran first and mutated `df_map`
     (`:478`, `:488-489`).
  3. `:531-534` `signal_cuts` =
     `minv > 1.08 && minv < 2.9 && pair_pt > 8 && FiducialGapCutExpr(m1.charge*m1.eta)
      && FiducialGapCutExpr(m2.charge*m2.eta)`, applied at `:549` (OS) and `:720` (SS).
  4. OS = `muon_pair_tree_sign2`.
  **There is NO dR cut** (`dr` appears only as a histogram axis and as input to eps_dR / eps_reco)
  and **no `q*eta < 2.2`** — replaced by the gap windows on 2026-08-17.
  `ParamsSet::single_mu_fiducial_gap_cuts` = {(-1.20,-1.05), (-0.10,+0.06), (+2.30,+2.40)}
  (`ParamsSet.h:489`), rejected CLOSED and compared in float (`ParamsSet.h:294-303`).

  **=> ANSWER: yes — the pp24 single-b signal region is exactly OS + Tight +
  `1.08 < m < 2.9 GeV` + `pair pT > 8 GeV` + the fiducial gap cut on both muons. Nothing else.**

- 2026-08-25 Step 1c — **F3 RESOLVED, and NOT in the direction the user expected. The Tight WP is
  NOT applied at the ntuple stage.** `requireTight = false` is the default and the pp24 run script
  does not override it, so the pair trees are **Medium**-quality with Tight recorded as a branch.
  Tight is applied **only inside `FillHistogramsCrossx`** (`PP.cxx:543-546`). The driver
  `RDFBasedHistFillingData::FillHistograms` (`Data.cxx:146-168`) calls
  `FillHistogramsGeneric()` at `:151-153` and only then `FillHistogramsCrossx()` at `:155-158`,
  so **every generic histogram is a MEDIUM-WP yield corrected by TIGHT eps_trig and TIGHT
  eps_reco** — a genuine mismatched correction, exactly the case
  `tight_wp_default_change.md` §4 forbids ("do NOT keep applying the Medium reco-eff to a
  Tight-selected spectrum" — here it is the mirror image, a Tight efficiency on a Medium
  selection). This is a REAL BUG in the generic family and is fixed by hoisting the
  `pair_pass_tight` filter above the generic booking (Step 5), not by touching the ntuple stage.

- 2026-08-25 Step 1d (scout: the generic booking and the jacobian defect, from code).
  Generic 1D list (`Data.cxx:372-378`), booked for each of `_ss`, `_op` with
  `generic_weight_col = w_reco_trig` (`PP.cxx:486`, `:491`):
  `h_Dphi`, `h_Dphi_zoomin`, `h_DR`, `h_DR_zoomin`, `h_Deta_zoomin`, `h_minv_zoomin`,
  plus a `_wgapcut` twin of each and `h_DR{,_zoomin}_<sign>_jacobian_corrected` (+ `_wgapcut`).
  Binnings from `var1D_pp.json` (`PP.cxx:9`): Dphi 64x(-pi,pi); Dphi_zoomin 20x(-0.8,0.8);
  DR 40x(0,5.75); DR_zoomin 20x(0,0.8); Deta_zoomin 20x(-0.8,0.8); minv_zoomin 40x(0,3.0).
  **No 2D or 3D generic histograms exist** — `df_filter_to_var2D/3D_list_map` is populated only
  for the trigger-efficiency filters, never for `_ss`/`_op`.
  The `_wgapcut` twins use `MuPairPassGapCut` (`Data.cxx:170-184`, `eta_gap_cut1 = 0.135` plus
  `charge_eta_gap_cuts`) — **a DIFFERENT and wider object than the analysis fiducial cut**, and
  since 2026-08-18 the comparison deliberately reads the PLAIN histograms, so the twins are dead
  weight (`output_gapcut_hists`).

  **F1 mechanism, confirmed at source.** `BaseClass.cxx:463`'s `(suffix, df, weight_col, vars)`
  overload takes the weight column VERBATIM; only the other overload (`BaseClass.cxx:425-461`)
  consults `weight_specifier_to_column_map` (`:447-452`). `Data.cxx:207-214` takes the former with
  `weight_col = w_reco_trig` and `suffix = category + "_jacobian_corrected"`, so the suffix is a
  NAME ONLY. Moreover `weight_specifier_to_column_map` is **never populated on the data side at
  all** (it is filled only in `RDFBasedHistFillingPythia.cxx:9`, `...Powheg.cxx:10`,
  `...PowhegFullsimSingleMuon.cxx:100`, and `...PowhegTruth.cxx:6`, the last being the only real
  jacobian column in the codebase, `weight_norm_over_truth_dr`). So the fix is TWO parts: define a
  `w_reco_trig_over_dr` column on the data side AND route the `_jacobian_corrected` booking
  through it.

- 2026-08-25 Step 2 DONE (shared binning registration, orchestrator-owned).
  * `ParamsSet.h`: added the FINE crossx pair-eta axis as a single source —
    `N_PAIR_ETA_CROSSX_BINS = 44`, `PAIR_ETA_CROSSX_MIN/MAX = -+2.4`, and the generated edge
    vector `pair_eta_crossx_bins` (45 edges). VERIFIED identical to the literal axis it
    replaces: max |edge difference| vs `TH1D(44,-2.4,2.4)` = **8.9e-16**, i.e. double rounding
    only. It is documented there as a DIFFERENT object from the 9 coarse physics panels
    `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`.
  * `RDFBasedHistFillingBaseClass.cxx::BuildHistBinningMapBaseCommon`: registered
    `hist_binning_map["pair_eta_crossx"]` plus the six 1D comparison axes
    `dr_bins_1d` (40x[0,5.75]), `dr_zoomin_bins_1d` (20x[0,0.8]),
    `dphi_bins_1d` (64x[-pi,pi]), `dphi_zoomin_bins_1d` (20x[-0.8,0.8]),
    `deta_zoomin_bins_1d` (20x[-0.8,0.8]), `minv_zoomin_bins_1d` (40x[0,3]).
    **Values UNCHANGED** — they are the pre-existing `var1D_pp.json` numbers, moved to one
    place so the data histogram and the MC partner of the SAME observable cannot be booked on
    different axes. That had been happening: fullsim `truth_dr_zoomin` 20x[0,1] against the
    data's 20x[0,0.8], `dphi/deta_zoomin` +-1.0 against +-0.8, `pair_eta` 48 bins against 44 —
    which is why no ratio pad was possible and every overlay was a density-only comparison.

- 2026-08-25 Step 8 (partial) — **six orphaned plot subdirectories deleted** (78 PNGs):
  `dphi/` (7), `hist_2D/` (25), `hist_2D_interval_projections/` (11, incl. `plots_with_4_to_8_GeV/`),
  `noDPHI_div/` (12), `powheg_data_compr/` (17), `with_gapcut/` (7).
  **KEPT pending an explicit user OK:** `backup/` (12 PNGs, 2023-03-07) and
  `backup_no_effcy_corr_w_mu4_mu4noL1/` (32 PNGs, 2026-06-10) — they meet the "no code can
  reproduce them" test only because they are DELIBERATE hand-made backups, which is a different
  thing from a stale output. The 12 top-level PNGs are kept until the new `signal/` and
  `generic/` sets exist, then removed.

- 2026-08-25 Step 4 DONE (MC-side RDF, subagent `_sub_mcRDF_1.md`).
  `var1D_pythia_fullsim.json` 24 -> 30 `variables1D`, six APPENDED, **no existing entry touched**:
  `truth_dr_zoomin_ppbin`/`truth_dphi_zoomin_ppbin`/`truth_deta_zoomin_ppbin`/
  `truth_minv_zoomin_ppbin` (binnings `d*_zoomin_bins_1d`), `truth_pair_eta_crossx`
  (`pair_eta_crossx`), `truth_pair_pt_log_150` (`pT_bins_150`).
  `RDFBasedHistFillingPythiaFullsim.cxx`: `BuildFilterToVarListMapExtra()` appends only (the
  reco-eff / det-response lists are untouched); `CreateBaseRDFsPythiaFullsimExtra()` hoists the
  existing `gap_truth`/`gap_reco` strings out of the category loop (pure refactor, same
  `ParamsSet::FiducialGapCutExpr` calls) and adds
  `df_{op,ss}_gapcut_truth_weighted = df_{op,ss}_weighted.Filter(gap_truth)`.
  **31 new keys**, weight column `weight` throughout. Signal family: the six 1D under each of
  `_single_b_pass_signal_truth`, `_op_pass_signal_truth`, `_ss_pass_signal_truth`, plus the 2D
  `h_truth_pair_eta_crossx_vs_truth_pair_pt_log_150_single_b_pass_signal_truth`
  (x = `pT_bins_150`, y = `pair_eta_crossx`) — the fullsim replacement for the truth sample's
  `h2d_sig_accept_num_pt_150_eta`. Generic family: the six 1D under `_op_gapcut_truth` /
  `_ss_gapcut_truth`.

- 2026-08-25 **PRE-EXISTING BREAKAGE FOUND AND FIXED AT SOURCE (orchestrator).**
  `RDFBasedHistFillingPythia.h` did not compile: line ~134 default-initialises
  `dr_bins_edges_for_reco_effcy` from `CommonEffcyConfig{}` but the header did not include
  `CommonEffcyConfig.h` -> "use of undeclared identifier 'CommonEffcyConfig'". Broken since the
  2026-08-18 10:13 edit that introduced that initialiser (the last `.so` predates it, 09:56),
  i.e. `RDFBasedHistFillingPythiaFullsim`, `...Truth` and `...FullsimOverlay` had ALL been
  uncompilable for a week. Fixed by adding the include to the header, with a comment recording
  why it must be there and not only in the .cxx.

- 2026-08-25 Step 3 DONE (data-side RDF, subagent `_sub_dataRDF_1.md`).
  a. **Tight WP hoisted (D3).** New `RDFBasedHistFillingPP::ApplyMuonWorkingPointFilter()`
     (`PP.cxx:451-487`), guarded by `wp_filter_applied` (`Data.h:257-261`), called first thing in
     the `!trigger_effcy_calc` branch of `FillHistogramsGeneric` (`PP.cxx:501`) and again at
     `PP.cxx:600` where the old inline block was. Exactly-once was proved over all five driver
     paths (generic on/off, `trigger_effcy_calc` true/false, `isTight=false`); `FillHistogramsCrossx`
     is reachable only when `trigger_effcy_calc == false`, which is precisely the branch that
     applies the WP, so the two callers cannot disagree. The Medium escape hatch is untouched.
  b. **Signal-region 1D set (§3a).** `PP.cxx:760-807` (OS) and `:849` (SS), filled from the SAME
     nodes and weight as the crossx 2D/3D (`df_single_b_crossx_weighted` / `df_ss_weighted`,
     `crossx_weight_trig_corr`), NOT width-scaled at fill time (matching
     `h1d_crossx_minv_0_4_op_dsigma`):
     `h1d_crossx_{DR_zoomin,Dphi_zoomin,Deta_zoomin,minv_zoomin,pair_eta,pair_pt_150}_w_signal_cuts_{op,ss}_dsigma`.
  c. **F1 fixed.** `PP.cxx:536-538` defines `w_reco_trig_over_dr = dr > 0 ? (w_reco*w_trig)/dr : 0`
     (form copied from `RDFBasedHistFillingPythiaTruth.cxx:96`); new
     `RDFBasedHistFillingData::BookJacobianCorrectedGeneric()` (`Data.cxx:186-220`) routes the
     `_jacobian_corrected` booking through it and THROWS if the jacobian column equals
     `generic_weight_col`, so the histogram can never again be a silent clone. Confirmed on 200k
     real entries: jacobian 585269 vs plain 161157 (previously bit-identical).
  d. **Task D.** The `44, -2.4, 2.4` literal is retired at **five** fixed-bin sites
     (`PP.cxx:760,763,772,793,818` -> `ParamsSet::N_PAIR_ETA_CROSSX_BINS/PAIR_ETA_CROSSX_MIN/MAX`,
     keeping them FIXED-bin so the on-disk representation is unchanged) and **three** edge-array
     sites (`:805,808,821` -> `pms.pair_eta_crossx_bins.data()`); both `make_unif_edges(44,...)`
     calls are gone.
  e. `var1D_pp.json`: the six 1D axes now request their binning BY NAME. Values unchanged; those
     six data histograms become variable- instead of fixed-binned in the output (edges identical
     to <= 1 ulp), which is what makes them bit-comparable with the MC partners.

- 2026-08-25 Step 6a DONE — **fullsim RDF rerun** (event loop 40.3 s wall / 204.7 s CPU, 8 slots).
  414 -> **445 keys**; 31 new, 0 lost. **Every pre-existing histogram's integral is unchanged**
  to 1e-9 relative (67 keys differ only in `GetEntries()` at the 1e-13 level — multithreaded
  summation order in the `_divided` graphs), so `pair_reco_eff_pp24_full.root` and the applied
  eps_reco are untouched and were deliberately NOT rebuilt.
  Signal-region MC integrals (nb; x1000 -> pb): single-b **7.84359**, all-OS 9.34088,
  all-SS 0.356349; the new 2D integrates to 7.84357, matching its 1D projection.
  So the signal-region comparison will read **MC 7844 pb vs data 5784 pb, MC/data = 1.356** —
  against 1.63 with the isospin-averaged truth sample, i.e. D1 removes ~17 % of the discrepancy
  and it was an artefact of comparing a Pb-isospin-averaged NN cross-section with pp data.

- 2026-08-25 Step 6b DONE — **pp24 data RDF crossx rerun** (event loop 8.6 s wall / 54.7 s CPU).
  47 -> **59 keys**; 12 new, 0 lost.
  * The **12 new signal-region 1D histograms** integrate exactly as predicted:
    OS `pair_pt_150` = **5784.0614** (bit-equal to `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts`),
    OS `pair_eta` / `DR_zoomin` / `Dphi_zoomin` / `Deta_zoomin` / `minv_zoomin` = **5784.0775**
    (= 5784.0614 + the 0.0161 that the pT_bins_150 axis pushes into overflow above 150 GeV);
    SS 502.6327 / 502.6281. That every non-pT observable lands on the SAME number is itself a
    check: no signal-region pair falls outside a zoom-in axis, as it must be, since
    m < 2.9 GeV with pair pT > 8 forces dR ~< 0.725 and hence |dphi|, |deta| < 0.725 < 0.8.
  * **The crossx spectra are UNCHANGED** (`h2d_crossx_*_w_signal_cuts` identical to 1e-6), which
    proves two things at once: the hoisted WP filter is applied exactly once (the crossx was
    already Tight), and the CURRENT `single_mu_effcy_pT_fit.root` (2026-08-24 17:14, F5) gives
    the same turn-on as the one the 2026-08-18 run used — the Aug-24 refit did not move pp24.
  * **The generic histograms drop by ~12 %** (`h_DR_op` 1.2727e7 -> 1.1233e7, ratio 0.883;
    every generic key moves by 0.875-0.895). That is the D3 fix: they are now Tight-WP yields
    instead of Medium-WP yields corrected by Tight efficiencies.
  * **F1 fixed and visible:** `h_DR_zoomin_op_jacobian_corrected` / `h_DR_zoomin_op` =
    9.26115e6 / 2.72855e6 = **3.394** (1/dR over [0, 0.8]) and
    `h_DR_op_jacobian_corrected` / `h_DR_op` = 1.28464e7 / 1.12334e7 = **1.144** over the full
    range, where the away-side dR ~ 3 pulls the weight below 1 — previously the two were
    bit-identical. (Corrected 2026-08-25 from the 3.00 / 1.01 first written here, which were an
    arithmetic slip; caught by the code review.)
  * Evaluator census printed by the run (H2's PrintStats): eps_dR 46.8 % on the exponential fit,
    0.018 % polynomial, 0.0012 % raw-bin placeholder, 53.2 % outside the cell grid (eps_dR = 1);
    eps_reco 87.3 % from the 3D map, 12.7 % from the dR-integrated fallback, 0 % inclusive,
    0 floored — with 52.1 % / 22.4 % / 73.5 % of evaluations clamped in pair pT / pair eta / dR.
    Those clamp fractions are the GENERIC domain (F4/O2) and are exactly what `generic/README.md`
    must warn about.

- 2026-08-25 Steps 5+7 DONE (plotting, subagent `_sub_plotting_1.md`).
  * `PlotMCDataComprBaseClass.{h,c}` rewritten: `output_subdir` + `OutputPath()` (mkdir -p)
    replaces the hardcoded save paths; the dead `pythia_with_data_resonance_cuts_dir` is gone;
    Pythia repointed to the fullsim full sample with the isospin reason in a comment;
    **`TRUTH_COMBINE_RENORM` deleted** so `norm_factor[pythia] = 1e3` (nb->pb) alone; a
    `data_already_dsigma` flag zeroes the 1/L factor for the signal family; ONE
    `McDataObservable` table is the only place a histogram name is spelled; `AssertSameAxis1D()`
    throws on any nbins/edge mismatch at relative tolerance 1e-9.
  * `plot_mc_data_compr_signal.cxx` (NEW) -> `signal/`, 5 panels, log y, per-pad MC selection
    named in the legend entry (OS pad = single-b, SS pad = inclusive SS).
  * `plot_mc_data_pair_pt_in_eta.cxx` repointed to the fullsim 2D, output -> `signal/`.
  * `plot_mc_data_compr.cxx` is now GENERIC-only, 8 panels, every `_unity` path removed (item 2).
  * **D6 verified bin by bin, not just on integrals.** The POWHEG flavour axis IS a partition:
    mechanism-sum vs flavour-sum agree to a max per-bin relative difference of 6.4e-15 (bb SS),
    1.4e-14 (bb OS), 8.1e-16 (cc SS), 2.6e-14 (cc OS). The old sum-of-ten was 28420.73 against
    the true 9864.32 (bb OS) = **2.881x**; 2.000x in the cc file, whose origin axis is empty.
    Recorded oddity: in the cc file the histogram NAMED `bb` holds the whole sample and `cc` is
    empty — a producer misnomer, harmless because the three flavour histograms are summed.
  * Two more latent compile breakages found and fixed at source by the orchestrator:
    `plotting_codes/helper_functions.c` used `std::cout` with no `<iostream>` (it only ever
    compiled when an includer happened to pull it in first), and see the
    `RDFBasedHistFillingPythia.h` entry above.
  * `pipelines/pipeline_pp_crossx.sh` Stage 8 now runs **three** macros instead of one;
    `plot_mc_data_pair_pt_in_eta.cxx` was manual-only and had therefore been going stale.

- 2026-08-25 Step 4b (orchestrator) — the generic family's full-range `DR` and `Dphi` panels had
  NO Pythia curve, because `McVsDataVar1Ds()` carries only the zoom-ins. Added
  `truth_dr_ppbin` (`dr_bins_1d`, 40x[0,5.75]) and `truth_dphi_ppbin` (`dphi_bins_1d`, 64x[-pi,pi])
  to `var1D_pythia_fullsim.json` and a `McVsDataGenericOnlyVar1Ds()` list attached to the generic
  filters ONLY — inside the signal region dR is kinematically capped at ~0.725 so a full-range
  view there would add nothing, whereas the generic family is precisely where the away-side peak
  at dR ~ pi lives. Fullsim RDF rerun a second time to produce them.

- 2026-08-25 Step 7 DONE — **both plot families produced.**
  `plots/mc_data_compr/signal/` (7 PNGs): `DR_zoomin`, `Deta_zoomin`, `Dphi_zoomin`,
  `minv_zoomin`, `pair_eta`, `pair_pt`, `pair_pt_in_eta_subplots`.
  `plots/mc_data_compr/generic/` (8 PNGs + `README.md`): `DR`, `DR_jacobian_corrected`,
  `DR_zoomin`, `DR_zoomin_jacobian_corrected`, `Dphi`, `Dphi_zoomin`, `Deta_zoomin`,
  `minv_zoomin`. The 12 superseded top-level PNGs were deleted.
  **Headline: data 5784.06 pb vs MC 7843.57 pb, MC/data = 1.356** (was 1.63 against the
  isospin-averaged truth sample).
  Two visual confirmations that the cut matching actually took:
  * `signal/pair_eta` shows the three gap-cut dips at eta^pair ~ 0, +-1.1, +-2.3 in **BOTH** the
    data and the MC, at the same positions — previously the MC had no gap cut at all and the
    dips were a data-only artefact of the comparison.
  * `signal/DR_zoomin` falls to exactly zero above 0.72 on both sides (data bin contents
    84.22 -> 28.50 -> 2.15 -> 0; MC 0.0693 -> 0.0224 -> 0.00136 -> 0), independently reproducing
    the kinematic bound dR ~< 2m/pT = 0.725 for m < 2.9 GeV at pair pT > 8 GeV.
  `generic/README.md` written per §4 constraint 5 / D4, with the clamp census refreshed to the
  numbers this run actually printed.
  `RDFBasedHistFillingPbPb.cxx` recompiled clean against the changed `RDFBasedHistFillingData.h`
  (additive changes only); PbPb outputs untouched and not rerun.

- 2026-08-25 Step 9 DONE (POWHEG, subagent `_sub_powheg_1.md`) — **additive, per D8.**
  `RDFBasedHistFillingPowhegTruth.cxx` gains a file-local anonymous namespace
  (`McVsDataGenericVar1Ds` / `McVsDataGenericFlavourSplit` / `McVsDataGenericSigns` /
  `kGapCutTag = "_gapcut_truth"` / `McVsDataGenericFilters`), a
  `df_{ss,op}_gapcut_truth_weighted = df_{ss,op}_weighted.Filter(gap_truth)` built from
  `ParamsSet::FiducialGapCutExpr("m1.truth_charge * m1.truth_eta")` and its m2 twin, and an
  inline generic-family fill block. **`FillHistogramsSignalAcceptance` and its `q*eta < 2.2`
  (now `:224-226`) are BYTE-UNCHANGED**, so the NLO template-fit input is untouched, as D8
  requires. `var1D_powheg_truth.json` append-only, 24 -> 30 (no `truth_` prefix — that is this
  file's convention, every variable in it is already truth).
  **60 new keys**: `h_<VAR>_{op,ss}_gapcut_truth` plus
  `..._flavor_binned_{single_b,both_from_b,both_from_c,others}`, VAR in
  {`dr_ppbin`, `dr_zoomin_ppbin`, `dphi_ppbin`, `dphi_zoomin_ppbin`, `deta_zoomin_ppbin`,
  `minv_zoomin_ppbin`}. The flavour partition sums to the inclusive to <= 4e-12 in all 12 cases,
  so D6's 2.881x over-count is now structurally impossible. All six axes verified EDGE-BY-EDGE
  equal to both the pp24 data and the Pythia fullsim (12/12 exact).
  Gap-cut cost, sign-independent as it must be: OS 14122.18 -> 11778.27 (**0.8340**),
  SS 5764.19 -> 4793.22 (**0.8316**) — the eps_acc-sized ~16.6 % offset.

  **Three findings that change what the plots may claim:**
  * **`weight_norm` is a genuine cross-section in pb, NOT nb** — `weight_norm = weight/N_gen`
    with `weight = EventWeights[0] * filter_effcy` (`PowhegAlgCoreT.c:209-212`,
    eps_filt = 0.003 bb / 0.001108 cc). Mean `EventWeights[0]` = 1.1e7; read as nb that would be
    an 11 mb b-bbar cross-section, 55x the total sigma_bb at 5.02 TeV — impossible; as pb it is
    11 ub, normal. **So the plotter's POWHEG normalization factor of 1 is CORRECT and POWHEG must
    NOT receive Pythia's nb->pb x1000.** Cross-check: POWHEG OS 14122 pb vs pp24 data OS 28055 pb.
  * **The POWHEG curve is bb-ONLY.** `muon_pairs_powheg_cc_truth.root` does not exist (only
    `_part{1,3,4,5}`; parts 2 and 6 were never produced) and the class merely prints
    "Skipping missing input file", so `others` is 0 everywhere. Hadding the incomplete cc parts
    would change the normalization and the flavour mix — a physics decision, deliberately NOT
    taken. The legend must say bb, not bb+cc.
  * **The on-disk POWHEG signal-acceptance histograms were STALE**, and the rerun refreshed them:
    4 of 400 pre-existing keys moved (`h2d_sig_accept_num_pt{,_150}_eta` + their two ratios).
    Not caused by this change — the file was dated 2026-04-13 and predated commit `7b74dfc`
    (2026-06-22), which deleted `&& truth_dr > 0.05` from `signal_cuts`. Proven by direct count:
    582146 signal pairs, 747 with `truth_dr <= 0.05`, and 582146 - 747 = 581399 = the OLD entry
    count exactly. Denominators identical to 1.4e-14. The other 396 are bit-identical.
  Also corrected in passing: the RDF stage writes ONE combined file
  `powheg_full_sample/histograms_powheg_truth.root`; the `{bb,cc}_evgen_truth_full_sample/`
  files the plotter still opened are 2025 legacy and are neither read nor written by the
  current producer.

- 2026-08-25 **BLOCKING BINNING BUG FOUND BY /review-plot AND FIXED BY USER DECISION:
  the 9 pair-eta panels OVERLAPPED and were MISLABELLED.**
  `plot_mc_data_pair_pt_in_eta.cxx:122-124` (and, identically, `SingleBCrossxPlotterBase.cxx`
  `:192-193,329-330` and `SignalAcceptancePlotter.cxx:117-118` — so the pp24 crossx plot set and
  the signal-acceptance plots too) project the 9 coarse panels out of the fine pair-eta axis with
  `FindBin(lo+1e-6) .. FindBin(hi-1e-6)`. With **44** uniform bins on [-2.4,2.4] the width is
  0.109090..., and **none** of the 8 internal panel boundaries (+-0.5, +-1.0, +-1.5, +-2.0) is a
  bin edge. Measured on the real histogram: bins **4, 9, 13, 18, 27, 32, 36, 41** each land in
  TWO adjacent panels, the 9 panels sum to **6942.93 pb** against the true total **5784.06 pb**
  (**+20.04 %**), and the panel labelled `[-1.0,-0.5]` actually draws `[-1.0909,-0.4364]` — a
  label/binning disagreement, i.e. exactly `.claude/CLAUDE.md` Binnings rule 5.
  **USER DECISION 2026-08-25: move the canonical fine pair-eta axis 44 -> 48 bins.** Width
  becomes exactly 0.1 and every coarse boundary is an exact edge (verified: max distance to the
  nearest edge = 4.4e-16 over all 8), so the overlap and the mislabelling disappear at source
  rather than being patched per plot. `sigma_fid = 5784.06 pb` is UNCHANGED — only the per-panel
  split moves. Rerun blast radius executed here: pp24 data RDF + pp24 fullsim RDF + all replots.
  **Pb+Pb deliberately NOT moved** (`RDFBasedHistFillingPbPb.cxx:1139,1294` still retype
  `44,-2.4,2.4`): leaving its code and its on-disk histograms mutually consistent is safer than
  changing the code without rerunning. Pb+Pb must adopt the constant AND rerun in one step.

- 2026-08-25 Step 10 — **48-bin rerun executed and the overlap VERIFIED GONE.**
  Fullsim RDF + pp24 data RDF rerun on `N_PAIR_ETA_CROSSX_BINS = 48`. Measured on the produced
  files, summing the 9 `pair_eta_proj_ranges_coarse_incl_gap` panels with the SAME
  `FindBin(lo+1e-6)..FindBin(hi-1e-6)` projection that produced the +20.04 % error:

  | | axis | total | panel sum | difference | bins in 2 panels | bins in none |
  |---|---|---|---|---|---|---|
  | data | 48 x 0.1000 | 5784.061 pb | 5784.061 pb | +1.8e-12 | 0 | 0 |
  | MC   | 48 x 0.1000 | 7843.570 pb | 7843.570 pb |  0.0     | 0 | 0 |

  `sigma_fid` and MC/data are unchanged (5784.06 pb, 7843.57 pb, **1.3561**), exactly as the
  decision predicted — the total never depended on the panel split. The 1D cross-checks also
  hold: `h1d_crossx_pair_pt_150_..._op` = 5784.0614 (bit-equal to the 2D) and
  `h1d_crossx_pair_eta_..._op` = 5784.0775 (the 0.016 pb difference is the pT > 150 GeV tail the
  `pT_bins_150` axis sends to overflow).
  NOTE for whoever next runs them: `plots/single_b_analysis/pp24{,_pt_150}/` and the
  signal-acceptance plots carried the SAME +20 % panel overlap and are fixed by this change too,
  but they have not been regenerated in this task — they are stale until the crossx plotting
  stage is rerun.

- 2026-08-25 Step 10b — **the pp24 crossx plot set regenerated on the fixed axis** (blast radius
  of the 44 -> 48 change, per the Autonomy Contract's "regenerate every result the fix affects"):
  `plots/single_b_analysis/pp24/` (5 PNGs), `plots/single_b_analysis/pp24_pt_150/` (3), and
  `plots/sanity_check_crossx/PP_2024_pair_pt_in_eta_subplots.png`. Their 9 pair-eta panels
  carried the SAME +20 % overlap and the same wrong labels via
  `SingleBCrossxPlotterBase.cxx:192-193,329-330`; both are now correct without any change to that
  plotter — the axis fix reaches it.
  `plots/sanity_check_crossx/PbPb_10-20%_pair_pt_in_eta_subplots.png` was re-rendered in the same
  run but is UNCHANGED in content: Pb+Pb code and Pb+Pb histograms are both still on the 44-bin
  axis, so **the Pb+Pb panels still carry the overlap**. That is the deliberate consequence of not
  moving Pb+Pb (Remaining Work).
  Docs updated at source: `docs/muon_wp_registry.md` §3 (it still claimed the data `isTight` was
  "declared but UNUSED ... no Filter", stale since `tight_wp_default_change.md` S2.3, and it is
  precisely the place this task changed) and `docs/signal_selection_change_impact.md` (a new
  header note giving the NARROWER rerun set for a change to this axis, and the explicit warning
  that Pb+Pb must switch the constant AND rerun in the same step).

## Results & Observations

(to be filled)

- 2026-08-25 Step 11 DONE (ratio pads + the remaining /review-plot fixes, subagent
  `_sub_plotfix_1.md`). Two new namespace-scope headers so the non-derived
  `plot_mc_data_pair_pt_in_eta.cxx` can share them: `McDataComprConfig.h` (WP enum + every input
  path builder, one place) and `McDataComprRatio.h` (pad split, ratio build, auto range,
  log-label policy, unity line, legend top-guard).
  * **Ratio pads on every 1D panel of both families and on both pair-pT views.** On the OS signal
    pads there are TWO ratio markers, single-b (black) and all-OS (blue), both over the data:
    their vertical separation IS the MC's own non-single-b OS content, which is the scale on
    which the data's unsubtracted background should be judged. The 9-panel figure stayed legible
    (cells 450x350 -> 450x450) and got ratio pads on one common range.
  * Jacobian y title now `(1/#DeltaR) d#sigma/d#DeltaR [pb]`; log-label policy generalised
    (`SetMoreLogLabels` under 3 decades, `SetNoExponent` under 1.5, a stricter 1.5-decade cut in
    the one-third-height ratio pad); legend headroom x6 plus a post-paint 0.88 NDC top guard, so
    no legend touches a frame line or a marker; `#DeltaR` stray space gone;
    `d#sigma/dp_{T}^{pair}`; 9-panel y range tightened by percentile (1.4e-4, not ~1e-5).
  * **POWHEG repointed** to `histograms_powheg_truth.root`, INCLUSIVE
    `h_<var>_{op,ss}_gapcut_truth`; the two POWHEG `DataType` entries collapse to one
    (`s_nDtTypes` 4 -> 3), and the flavour summing plus `GetPowhegHist1D` are DELETED — D6's
    2.881x over-count is now impossible by construction, not by discipline. Verified: flavour
    sum 11778.3 vs inclusive 11778.3 (3.6e-12). It now passes through `AssertSameAxis1D` and
    into the ratio, and the legend reads "POWHEG b b-bar" because the sample is bb-only.
  * **WP config var added** (repo rule): macro argument `"tight"|"medium"`, default Tight. It
    switches the DATA file to the `_medium_wp` crossx output (selection + corrections together)
    and CANNOT switch the MC, every MC histogram here being a truth quantity for which a
    reconstruction WP does not exist. Registered in `docs/muon_wp_registry.md` §3b.
  * **Two ordering bugs caught only by LOOKING at the rendered PNGs**, both of which would have
    shipped silently: `MakeRatio` clones the numerator, so a ratio built after `HideXAxis()`
    inherited the suppressed axis and the generic ratio pads came out with NO x axis at all; and
    the data-only jacobian panels were given an empty ratio pad and kept ROOT's 10 % left margin
    in the un-split branch, which deleted the y title and chopped the leading digit off every y
    label. Ratios are now built before any main-pad styling in all three macros.
  * Reference numbers reproduced exactly after every change: data 5784.06 pb, MC 7843.57 pb,
    MC/data 1.35607, panel projection exact, no `[SKIP]` or `[WARN]` on any run.

- 2026-08-25 **NEW PHYSICS FINDING — POWHEG and Pythia do NOT carry the same low-mass selection
  in the generic family.** The plotter showed POWHEG's OS `minv` identically zero below 1.05 GeV,
  exactly like the data and unlike Pythia. Confirmed at source: `PowhegAlgCoreT.c:298-350` builds
  `resonance_tagged_muon_index_list` and drops any pair containing a tagged muon
  **UNCONDITIONALLY**, with no `turn_data_resonance_cuts_on` switch, whereas
  `PythiaAlgCoreT.c:910,1013` gates the same veto behind that flag and the Pythia file in use is
  the `_no_data_resonance_cuts` build. So in `generic/` the three curves carry THREE different
  low-mass selections: data vetoed, POWHEG vetoed, Pythia not. Inert inside the signal window;
  recorded in `generic/README.md` and in Remaining Work.

- 2026-08-25 Housekeeping: `plot_mc_data_2D_hists_and_1D_proj.cxx` DELETED (recoverable from git
  history). It could not run before this task (three independent blockers, Step 1a), both
  directories it wrote were deleted as orphaned, and after the POWHEG `DataType` collapse it no
  longer compiles either — leaving a non-building dead file in the tree is a trap, and porting it
  would mean rewriting physics into a file that can never execute.
  `docs/muon_wp_registry.md` gained §3b (the plot set's WP config var) and had its stale §3
  corrected.

### Review findings and their disposition (2026-08-25)

Both mandated reviews ran as read-only subagents. `/review-analysis-code` -> **FAIL** on one build
defect (now fixed), PASS on all five physics questions; `/review-plot` -> **PASS-WITH-COMMENTS**,
no CRITICAL.

| # | sev | finding | disposition |
|---|---|---|---|
| C1 | HIGH | The signal-family colour constants were `protected` members of `PlotMCDataComprBaseClass`, but `plot_mc_data_pair_pt_in_eta.cxx` is a free-function macro that never includes the class -> ACLiC fails, `root` exits 1, and Stage 8 runs under `set -Eeuo pipefail`, so **the generic family would never have been produced at all**. | FIXED: constants moved to namespace scope in a new `McDataComprColors.h` included by both producers. |
| C2 | LOW-MED | `wp_filter_applied` is an incomplete ledger — `PP.cxx:151` applies the same filter in the trigger-efficiency path without recording it. | ACCEPTED, documented: the two paths are mutually exclusive (`trigger_effcy_calc`), so the ledger cannot be wrong where it is read. |
| C3 | MED | The new MC generic booking was wrapped in the file's local `try`/`catch`-and-continue idiom — the silent-incomplete-output shape of `reference_root_swallows_rdf_exceptions`: the MC curve would vanish and the plotter would print a `[SKIP]` nobody reads. | FIXED: the block now throws, with `map_at_checked` naming the missing key. |
| C5 | LOW | The Progress Log's jacobian ratios (3.00 / 1.01) did not match the file. | FIXED: they are **3.394** and **1.144**; an arithmetic slip of mine, corrected above. |
| P1 | HIGH | The 9 pair-eta panels overlapped and were mislabelled (+20.04 % on the panel sum). | FIXED AT SOURCE by user decision: fine pair-eta axis 44 -> 48 bins. See the Step-10 entry. |
| P2 | HIGH | No ratio panel on any canvas, which R3 makes mandatory for a comparison figure — and identical axes (D2) existed precisely to make one possible. | FIXED by user decision: MC/data ratio pads on BOTH families. |
| P3 | HIGH | The two `*_jacobian_corrected` panels carried a y title byte-identical to the unweighted panel beside them, so the two were indistinguishable from the image. | FIXED: title now names the 1/dR weighting. |
| P4 | MED | POWHEG drawn on a different binning, never axis-checked, normalization assumed. | FIXED by Step 9: POWHEG is on the shared axes (12/12 edge-exact), enters the assertion, and its factor-1 normalization is now DERIVED (`weight_norm` is in pb) rather than assumed. Remaining: it is bb-only. |
| P5 | MED | `generic/Deta_zoomin` SS spans 0.82 decades on a log axis and therefore shows ZERO y-axis numeric labels. | FIXED: `SetMoreLogLabels()` / linear fallback under ~1 decade, applied as a general rule. |
| P6 | MED | `generic/README.md` accurate but incomplete. | FIXED: census refreshed to this run, plus the near-side MC/data ~ 1038 point, the resonance-veto zero regions, POWHEG, and the data-only jacobian panels. |
| P7 | MED | No `signal/README.md`, while the generic one points readers there for anything quantitative. | FIXED: written, with the fiducial / background-unsubtracted / non-unfolded and pure-truth-single-b caveats. |
| P8 | MED | The `minv_zoomin` axis (40 x [0,3]) has no edge on 1.08 or 2.9, so the first and last populated bins are 60 % / 66.7 % inside the window and `Scale(,"width")` shows an artificial drop. | DOCUMENTED in `signal/README.md`, not fixed: fixing it means another canonical binning change, which is a separate user decision. Carried in Remaining Work. |
| P9 | MED | `docs/muon_wp_registry.md` §3 stale in exactly the place this task changed; and this plot set exposes no Medium/Tight WP config var. | Registry FIXED at source; the WP config var added to the macros. |
| P10/P11 | LOW | Legend on the frame line / touching markers; `#Delta R` renders with a stray space; `pair_pt` y title does not name the pair pT; one global y range with dead space. | FIXED. |

## Remaining Work

1. **`backup/` and `backup_no_effcy_corr_w_mu4_mu4noL1/` were NOT deleted** — they are deliberate
   hand-made backups, not stale outputs, so they fail the "no code reproduces them" test for a
   different reason than the six directories that were removed. One word from the user removes them.
2. **Two pair-eta binnings still coexist in the pp24 output, both pre-existing and both outside
   this task's scope.** `RDFBasedHistFillingPP.cxx:699-700`
   (`h2d_crossx_minv_0_4_vs_pair_eta_{op,ss}_dsigma`, the low-mass template block) retypes
   `24, -2.4, 2.4`, against the 44-bin `pair_eta_crossx` everything else now shares. Someone must
   decide whether the template block should adopt the canonical axis — it is a template-fit input
   (`low_mass_dimuon_template_fit.md`), so changing it has a blast radius there.
3. **Pb+Pb still on the OLD 44-bin pair-eta axis** (`RDFBasedHistFillingPbPb.cxx:1139,1294`
   retype `44, -2.4, 2.4`). Its code and its on-disk histograms are mutually CONSISTENT, so
   nothing is silently wrong today — but its 9 pair-eta panels still carry the **+20 % overlap and
   the wrong labels** that pp24 has just been cured of (visible in
   `plots/sanity_check_crossx/PbPb_10-20%_pair_pt_in_eta_subplots.png`). Pb+Pb must switch to
   `ParamsSet::N_PAIR_ETA_CROSSX_BINS` **and** rerun its crossx hist filling in the SAME step;
   doing one without the other recreates exactly the silent mismatch this task removed.
   Deliberately not done here: Pb+Pb is out of scope, its signal region still differs
   (`q*eta < 2.2`), and R_AA must not be recomputed until it is brought over.
4. **THREE different low-mass selections coexist in the generic family.** The data OS tree is
   resonance-vetoed at the ntuple stage; **POWHEG is too** (`PowhegAlgCoreT.c:298-350`, applied
   UNCONDITIONALLY — no `turn_data_resonance_cuts_on` switch); **Pythia is not** (the file in use
   is `_no_data_resonance_cuts`). Fixing the Pythia side needs a full MC NTuple Condor rerun of a
   4.5 GB tree. Documented in `generic/README.md` instead. Inert inside the signal window.
   Worth confirming on the POWHEG producer side that the unconditional veto is intended.
5. **POWHEG: RESOLVED (Step 9), with one residue.** It now carries the same fiducial gap cut as
   the data and sits on the shared 1D axes (12/12 edge-exact), so it enters the axis assertion and
   the ratio like any other curve, and `FillHistogramsSignalAcceptance`'s `q*eta < 2.2` — the NLO
   template-fit input — is byte-unchanged. RESIDUE: the sample is **bb-only**
   (`muon_pairs_powheg_cc_truth.root` does not exist, parts 2 and 6 were never produced), so a
   cc contribution is simply absent. Producing the missing cc parts, or hadding the incomplete
   ones, changes the normalization and the flavour mix and is a physics decision left to the user.
6. **POWHEG appears on the generic `DR` panel only** — its Mar-2025 file holds
   `h_DR_sign*_<flavour>` and nothing else (no `DR_zoomin`, `Dphi`, `Deta_zoomin`, `minv_zoomin`,
   no 1/dR-weighted variant). A file limitation, not a code one; it also predates the current
   producer naming (`_origin_binned_*` / `_flavor_binned_*` in `RDFBasedHistFillingPowhegTruth.cxx`)
   and still uses the one-sided `q*eta < 2.2` (`:158`) rather than the gap cut, so its selection
   does NOT yet match the data's. Rerunning POWHEG truth with the gap cut is the way to close
   both; not done (generic family is explicitly low priority).
7. **The `minv_zoomin` axis has no edge on the mass-window boundaries** 1.08 / 2.9 (it is
   40 x [0,3]), so the first and last populated bins are only 60 % / 66.7 % inside the window and
   `Scale(,"width")` shows an artificial drop at both ends (142.97 and 144.96 against a ~250
   plateau). A plotting artefact, not physics; documented in `signal/README.md`. Fixing it means
   another canonical binning change -> a separate user decision.
8. **The two `*_jacobian_corrected` generic panels are data-only.** No 1/dR-weighted histogram
   exists under the MC generic filters, and POWHEG's file has none either — drawing an unweighted
   MC curve against a weighted data curve would be F1 in mirror image.
9. **A Medium-WP data crossx output does not exist**, so the plot set's new WP config var throws
   on `"medium"`. Produce it with `RDFBasedHistFillingData::isTight = false` when the WP
   systematic is wanted.
10. The fullsim has ~4-7x fewer events than the truth sample in the top two pT-hat slices, so the
   50-150 GeV pair-pT tail of the signal plots is statistically weaker than it was. Visible as
   error bars; quote it with the high-pT points.

## Latest Stage

**DONE 2026-08-25.** All six Autonomy-Contract Done items are met. Both plot families are on disk
with ratio pads; both mandated reviews ran and every finding is closed or explicitly carried into
Remaining Work; the six implementation scratch docs are merged into this doc and deleted; the work
is committed. Nothing is in flight.

Deliverables at their canonical paths:
* `plots/mc_data_compr/signal/` — 7 PNGs + `README.md`
* `plots/mc_data_compr/generic/` — 8 PNGs + `README.md`
* `plots/single_b_analysis/pp24{,_pt_150}/` and `plots/sanity_check_crossx/PP_2024_*` regenerated
  on the corrected pair-eta axis (blast radius of the 44 -> 48 change).