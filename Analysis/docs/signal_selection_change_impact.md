# Signal-Selection Change — Impact & Rerun Map

**Type:** Repo reference doc (durable). **Scope:** what must be recompiled,
rerun, and re-plotted whenever the **single-b dimuon signal-region selection**
changes — a cut added, removed, or its value changed. Applies to *any* such
change, not one specific edit; it is also the checklist for **systematic
uncertainty** variations of the selection (e.g. ΔR, minv, pair-pT, q·η bounds).

> **First executed run of this map:** the ΔR>0.05 removal (2026-06-22) — see
> `docs/tracking/remove_dr_cut_signal_selection.md` for a worked example of the
> full chain below (recompile 7 classes → rerun pp+pbpb crossx + pythia truth →
> replot crossx/R_AA/acceptance/cutflow). Key lesson logged there: a cut can live
> in MORE than one signal-region definition per file (PythiaTruth had both the
> acceptance `signal_cuts` AND the template-fit `kin_cuts`) — always grep, don't
> trust an enumerated line list.
>
> **Ground truth for the selection itself:** `docs/analysis_overview.md` §2
> (signal region) and the relevant tracking-doc Physics Procedure. This doc
> owns only the *engineering dependency graph*. Pipeline/stage details:
> `README.md` + `pipelines/`.

> **A binning change has its own, smaller blast radius — and one such change landed
> 2026-08-25.** `ParamsSet::N_PAIR_ETA_CROSSX_BINS` (the FINE pair-eta axis of the crossx 2D/3D
> views) moved **44 -> 48**, i.e. width 0.109090... -> exactly 0.1. This is not a selection
> change: no cut moved, and `sigma_fid` is unchanged. It was necessary because with 44 bins
> **none** of the 8 internal boundaries of the 9 coarse panels
> (`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`) was a bin edge, so every
> `FindBin(lo+eps)..FindBin(hi-eps)` panel projection shared a bin with its neighbour — the 9
> panels summed **+20.0 %** above the true total and the panel LABELS disagreed with the bins
> drawn by up to 0.091 in eta. Affected consumers: `SingleBCrossxPlotterBase`,
> `SignalAcceptancePlotter`, `plot_mc_data_pair_pt_in_eta.cxx`.
> Rerun set for a change to this axis (much narrower than a selection change): the hist-filling
> stages that BOOK the axis (pp24 data crossx, pp24 Pythia fullsim) plus everything in §5 that
> projects it. Corrections are NOT affected — `eps_reco` is celled on
> `ParamsSet::pair_pt_coarse_bins` x `pair_eta_proj_ranges_coarse_incl_gap` x
> `dr_bins_edges_for_reco_effcy`, none of which is this axis.
> **Pb+Pb ADOPTED it on 2026-09-08** (`docs/tracking/mu_pt45_gap125_pairpt9_adoption.md` D1):
> `RDFBasedHistFillingPbPb.cxx` no longer retypes `44, -2.4, 2.4` in either spelling and now
> reads `ParamsSet::N_PAIR_ETA_CROSSX_BINS`, together with the fiducial + pair-level gap cuts,
> in the SAME step -- which is the condition this note demanded. Its crossx hist filling MUST be
> rerun before any Pb+Pb panel plot or R_AA is quoted, or its code and its on-disk histograms
> disagree.
> Owning doc: `docs/tracking/mc_data_compr_signal_generic_split.md`.

---

## 0. The current signal region (reference)

> **[!] pp AND Pb+Pb SHARE A SIGNAL REGION AGAIN (since 2026-09-08).** Pb+Pb was migrated onto the
> detector-gap fiducial cut, the pair-level `|eta^pair| < 2.2` window and
> `ParamsSet::N_PAIR_ETA_CROSSX_BINS` (48) in one step, together with the muon-pT and pair-pT
> changes below. **R_AA is formable again ONCE BOTH SIDES HAVE BEEN REFILLED** -- the code is
> migrated, the histograms are not yet. Owning doc:
> `docs/tracking/mu_pt45_gap125_pairpt9_adoption.md` (D1). For the 2026-08-18..2026-09-08 period,
> when the two genuinely differed, see `docs/tracking/pp24_crossx_rerun_2026_08.md`.

**pp24 AND Pb+Pb (reco), current (cut set changed by the user 2026-09-08):**
```
ParamsSet::SignalMinvCutExpr()  (= minv > 1.08 && minv < 2.9; signal_minv_min/max, SINGLE SOURCE since 2026-09-17)
  && pair_pt > ParamsSet::signal_pair_pt_min   <-- 9 GeV, was 8
  && BOTH muons pass ParamsSet::PassSingleMuFiducialGap(eta, charge)
  && |pair_eta| < ParamsSet::pair_eta_fiducial_max          <-- NEW, pair level
```

Plus, applied UPSTREAM at the **NTuple stage** and therefore inherited by every consumer:
**muon reconstructed pT > 4.5 GeV** (was 4.0), and the truth-pT analogue in MC. That cut is
**mode-dependent in data**: the trigger-efficiency NTuple mode (`trigger_effcy_calc`)
deliberately keeps **4.0 GeV** so the tag-and-probe retains the turn-on rise -- `eps^nc` is a
per-muon efficiency and is only ever evaluated above 4.5 (decision D8). MC is 4.5 everywhere,
with no loose variant (D9). Because this one lives at the NTuple stage it BREAKS the §1
"NTuple processing unchanged" boundary -- see the carve-out there.
i.e. q*eta = charge*eta must NOT lie in any window of
`ParamsSet::single_mu_fiducial_gap_cuts` = {(-1.25,-1.05), (-0.10,+0.06), (2.20,2.40)},
rejected on CLOSED intervals. Combined with the ntuple-level |eta| <= 2.4 the surviving region is
exactly `[-2.4, 2.20)` minus the two interior windows -- precisely the region the single-muon
turn-on fits cover (`CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap`, top edge moved
2.30 -> 2.20 WITH the cut; the coupling is enforced by a startup throw in
`RDFBasedHistFillingData::SetIOPaths`). That coupling is not cosmetic: a muon outside the fitted
range has NO efficiency and `EvaluateSingleMuonEffcyPtFitted` throws. **Both the predicate and
the RDF/JIT string form are built in FLOAT** (`ParamsSet::FiducialGapCutExpr`) so the Filter and
the float-valued binning agree bit for bit -- a double-vs-float mismatch at the forward edge is
exactly what broke the first 2026-08-18 run.

The **PAIR-LEVEL** window `|eta^pair| < 2.2` (symmetric, strict, `ParamsSet::
PairFiducialEtaCutExpr`) is new on 2026-09-07 and travels with the single-muon windows
EVERYWHERE they are applied to both legs of a pair -- signal region, generic family, template
inputs, pp24 fullsim reco-eff (truth AND reco legs), POWHEG truth/fullsim signal families, and
MC trig-eff Steps 2/3/4. It is NOT applied to MC trig-eff Step 1 (a single-muon map) nor to the
data tag-and-probe probe leg (user decision: eps^nc is per-muon). Motivation and the full
site inventory: `docs/tracking/muon_gap_cuts_acceptance.md` F17.

**Binning moved with it:** `CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` outer bins
+-(2.0, 2.4) -> +-(2.0, 2.2), so the 9 coarse pair-eta panels exactly tile the surviving region
(user decision, 2026-09-07). Consumers: pp24 + fullsim crossx hist filling, the pair reco-eff
cells, the MC trig-eff pair-eta histograms, and every plotter that projects those panels.

> **[RESOLVED IN CODE 2026-09-08 — historical note, kept because the guard it describes is still
> live.]** As found on 2026-09-07, Pb+Pb booked a RETYPED `44, -2.4, 2.4` pair-eta axis
> (width 0.10909...) on which **2.2 is not a bin edge**, and had adopted neither the fiducial cut
> nor the pair-level one. With the panel top at 2.2 the panel projection's `FindBin(hi - 1e-6)`
> landed on bin 43 (upper edge **2.29091**) instead of bin 44, so every Pb+Pb pair with
> |eta^pair| in [2.29091, 2.4] silently vanished from the panels and from the pair-eta-integrated
> dsigma/deta -- pairs no Pb+Pb cut removes -- while the outer panel, labelled (2.0, 2.2),
> actually integrated 1.96364-2.29091.
> **Still live, and the reason this note is kept:** `Utilities/PairEtaPanelBins.h`
> (`PairEtaPanels::Bins` / `CheckAxisAligned`, used by all ten pair-eta panel projection sites and
> by the two MC-trig-eff closure plotters) **THROWS** on any panel edge that is not a bin edge of
> the histogram's eta axis, so a mismatch fails loudly instead of silently.
> **What changed:** decisions D1 and D7 (`docs/tracking/mu_pt45_gap125_pairpt9_adoption.md`) moved
> Pb+Pb AND the Pythia/POWHEG truth producers onto `ParamsSet::N_PAIR_ETA_CROSSX_BINS` (48 over
> [-2.4, 2.4], width 0.1, on which every coarse panel edge IS a bin edge) and onto the fiducial +
> pair-level cuts. Pb+Pb is therefore **inside** the rerun set, not outside it.
> **What is still stale is the HISTOGRAMS ON DISK, not the code:** they stay on the retired 44-bin
> axis until the Pb+Pb crossx refill completes, so the guard will still throw on them and the
> Pb+Pb panels must not be plotted until then. That is why `pipeline_pp_crossx.sh` keeps
> `INCLUDE_PBPB_SANITY=false` by default.
> Found by `/review-analysis-code` 2026-09-07 (CRITICAL); owning doc
> `docs/tracking/muon_gap_cuts_acceptance.md` F17.

**Pb+Pb crossx and ALL truth analogs (Pythia/Powheg): MIGRATED 2026-09-08** onto exactly the form
above (decisions D1 and D7). The retired one-sided `m.charge*m.eta < 2.2` no longer appears in any
live selection in `RDFBasedHistFilling/`, `NTupleProcessingCode/` or `Utilities/`.
Truth analog: same with `truth_*` variables + `from_same_b` and `truth_pt > 4.5`.

**No dR cut anywhere** (removed 2026-06-22).

> **Note on q*eta:** the retired cut was one-sided per muon (`q*eta < 2.2`, no explicit lower
> bound; the floor is the muon |eta| ~ 2.5 detector edge). With the 2026-09-07 forward edge back
> at 2.20 the swap is yield-neutral THERE, and the cut's whole cost is the three gap windows:
> -6.93 % of pp / -9.03 % of PbPb reconstructed muons (Tight). **Since 2026-09-17 the gap cuts
> are ACCEPTANCE cuts carried by the pair reco efficiency** (`docs/tracking/pair_reco_eff_gap_acceptance.md`):
> they sit on the RECO leg of eps_reco only, the truth denominator carries the signal cuts
> alone, and the corrected pp24 cross-section refers to the truth region WITHOUT the single-muon
> q·η gap windows (muon p_T > 4.5 GeV, |η| < 2.4, mass window, pair p_T > 9 GeV) and with
> |η^pair| < 2.2 as the measured range — the truth fiducial edge of the eps_reco η^pair axis
> (`analysis_overview.md` §2). There
> is NO separate eps_acc; the historical values 0.8789 / 0.8765 are superseded and must not be
> applied. **Consequence for this map:** changing a gap window (single-muon or pair-level) is a
> change of eps_reco's numerator selection -> rerun the fullsim RDF hist filling, rebuild
> `pair_reco_eff_pp24_full.root`, refill the pp24 crossx -- in addition to the trigger-efficiency
> reruns already listed.
> The percentages quoted just above (-6.93 % / -9.03 %) were measured on the same superseded
> window and are stale for the same reason.

> **Note on the ΔR cut (motivation, for systematics):** `dr > 0.05` was added
> because the **data-based** dR-dependent trigger-efficiency inverse-weighting
> lacked statistics at small ΔR (the `minv > 1.08` cut drives the ΔR
> distribution toward zero at small ΔR). It induces a strong pair-pT–dependent
> signal acceptance (collimation: ΔR shrinks as pair pT rises — see the cutflow
> `…/single_b_analysis/pythia/pythia_sig_accept_above_60GeV*.png`). The dR
> trigger correction is being **moved to MC**; the correct (possibly absent)
> pair-pT–dependent ΔR cut can only be fixed once that MC exists.

---

## 1. What does NOT change (boundaries of the blast radius)

- **NTuple processing / `muon_pairs_*` trees** — UNCHANGED *by a signal-selection cut change*.
  Signal-region cuts are applied downstream in RDF hist-filling, not at tree creation, so a
  change to a cut listed in §0 does **not** require reprocessing ntuples or resubmitting Condor.

  > **[!] CARVE-OUT — this boundary is about SIGNAL-SELECTION cuts only, and it does NOT hold
  > for a change to the cuts that ARE applied at tree creation** (muon quality/WP, pT > 4.5,
  > |η| < 2.4, one-sided Δp/p, the impact-parameter cut, trigger matching, the resonance veto).
  > On **2026-09-08** the **muon pT threshold moved 4.0 -> 4.5 GeV** at this very stage
  > (`docs/tracking/mu_pt45_gap125_pairpt9_adoption.md`), in data AND MC, so **every** NTuple
  > stage had to be reprocessed -- pp24 both modes, Pb+Pb 2023/2024/2025/2026, the Pythia fullsim,
  > the HIJING overlay, and the Pythia/POWHEG truth productions. It is mode-dependent in data
  > (the trigger-efficiency mode stays at 4.0, decision D8). Note this also drags POWHEG
  > fullsim back into the blast radius, which the all-vertex change below had carved out:
  > that carve-out rested on the change being reco-only, and a truth-pT cut is not.
  >
  > On **2026-09-08** the pp impact-parameter cut became an **all-vertex, same-vertex** pair
  > requirement (`docs/tracking/pp24_all_vertex_pairs.md`), which is exactly such a change: the
  > pp24 trees themselves move, so **both** pp24 Condor modes (nominal `run_pp_24_nominal.sub`
  > and trigger-efficiency `run_pp_24.sub`) MUST be resubmitted, the ε^nc fits refit on the new
  > trees, and the pp-conditions **fullsim MC** NTuple stages rerun as well (the MC mirrors the
  > selection, so ε_reco moves with it). Pb+Pb is untouched by that change.
  >
  > The MC side is wider than ε_reco: the gate is the whole `FullSimSampleType::pp` production,
  > so the **detector response** and the **template-fit MC** move with it too (all four come out
  > of the same `pipeline_pythia_fullsim_pp.sh` Stage 5). The POWHEG **pp17** fullsim is also
  > switched on by its `!isFullsimOverlay` gate, but is deliberately NOT rerun — its only live
  > product is truth-level; see `docs/tracking/pp24_all_vertex_pairs.md`.
- **Trigger-efficiency derivation (P2 ε^nc fits, P3)** — UNCHANGED *by a pair-cut change*. The
  single-muon tag-and-probe ε^nc(pT, q·η) fits and their `_fine_q_eta_bin`
  inputs are independent of the *pair* signal-region cuts. `pipeline_pp_trig_eff.sh`,
  `pipeline_pbpb_trig_eff.sh`, `SingleMuEffcyPtTurnOnFitter` and all trig-eff
  plots stay valid. (The *applied* per-pair trigger weight inside crossx is
  recomputed automatically when crossx re-runs — no separate action.)

  > **[!] CARVE-OUT — this does NOT hold for a change to the SINGLE-MUON gap windows or to the
  > coarse q·η binning, and both changed on 2026-09-07.** The T&P **probe** leg applies
  > `ParamsSet::single_mu_fiducial_gap_cuts` (`PP.cxx:172`, `PbPb.cxx:340`), and the turn-on fits
  > are binned on `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap`, whose top edge is slaved
  > to the forward window and moved 2.3 → 2.2. So **P2/P3 MUST be refit** — it is the FIRST item
  > in the rerun order of `docs/tracking/muon_gap_cuts_acceptance.md` Latest Stage. Only a change
  > confined to *pair-level* cuts leaves the trigger efficiency alone.
- **Reco-eff PLACEHOLDER file** (`run2_reco_eff_placeholder.root`, Run 2 TF1,
  function of pT & q*eta only -- **no dR dependence**) -- still what **Pb+Pb** folds in as
  `w_reco`. **NOT pp24 any more:** since 2026-08-18 pp24 applies the genuine fullsim PAIR
  efficiency `pair_reco_eff_pp24_full.root` (`Utilities/PairRecoEffEvaluator.h`) plus the MC 2mu4
  dR-correlation correction (`Utilities/DrCorrectionCrossxEvaluator.h`), so for pp §3.C below is
  RESULT-AFFECTING, not bookkeeping. See §4.

---

## 2. [RECOMPILE] — code carrying the signal cut (ACLiC `.L file.cxx+`)

**Since 2026-09-17 the MASS WINDOW is no longer retyped anywhere live**: every site below reads
`ParamsSet::SignalMinvCutExpr("minv")` / `("truth_minv")` (and the cutflow macro builds its two
mass steps from `signal_minv_min/max`), exactly as the pair-pT threshold has read
`SignalPairPtCutExpr` since 2026-09-08. A mass-window change is therefore ONE edit in
`ParamsSet.h` + the recompile/rerun map below. **Extra blast radius of a mass-window change:**
the pp24 single-value pair trigger efficiency is MEASURED inside the window and applied above
74 GeV (`docs/tracking/pp24_trig_eff_hybrid_application.md` D3) — rerun
`pipelines/run_mc_trigeff_pair_eff.sh` (`PairTrigEff::Windows()` follows ParamsSet automatically;
`CheckSignalWindowMirror` throws until the file is regenerated) BEFORE the pp24 crossx refill.
Still retyped (dead / legacy code only): `SingleBAnalysis/*`, `OldMuonPairHistFillingCodePreRDF/*`,
`FillHistogramsCrossx_PP_clean.cxx`, `plot_reco_distr_singleb_vs_op_pp24.C` (STALE-gated).

Edit the cut in **all** of these (keep them in sync):

| File | Where | Feeds |
|---|---|---|
| `RDFBasedHistFilling/RDFBasedHistFillingPP.cxx` | :382 (crossx); :552 (ΔR-binned "no-minv" variant) | pp data crossx |
| `RDFBasedHistFilling/RDFBasedHistFillingPbPb.cxx` | :920 (crossx); :1195 (ΔR-binned variant) | PbPb data crossx |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaTruth.cxx` | :412 | truth signal acceptance |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegTruth.cxx` | :158 | Powheg truth acceptance |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsim.cxx` | `pass_signal_truth`/`pass_signal_reco` in `CreateBaseRDFsPythiaFullsimExtra` (was cited as :129/:131; :289/:291 as of 2026-09-03) | pp reco-eff / det-resp |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsimOverlay.cxx` | :69/:71 | PbPb (overlay) reco-eff |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegFullsim.cxx` | `pass_signal_truth`/`pass_signal_reco` in `CreateBaseRDFsPowhegFullsimExtra`, MIGRATED off the retired one-sided `q*eta < 2.2` on 2026-09-08 (D7) (was cited as :185/:186; :259/:260 as of 2026-09-03) | Powheg reco-eff (obsolete demo) |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegFullsim.cxx` | the `_single_b_pass_signal_truth_gapcut` filter in `CreateBaseRDFsPowhegFullsimExtra` (CURRENT cut set, from `ParamsSet`) | POWHEG FullSim pp17 curve of `plots/mc_data_compr/signal/pair_pt{,_in_eta_subplots}_mc_data_compr.png` — rerun the POWHEG fullsim RDF stage, then the two macros |
| `plotting_codes/single_b_analysis/plot_sig_accept_cutflow_above_60GeV.cxx` | :40 (`kCuts`) | cutflow diagnostic — **must mirror the cut list/order** |
| `Utilities/MCTrigEffPairSelection.h` | `SingleBSignalCutsReco()` | the **data-like mirror** of the PP `signal_cuts`, consumed by `FillMCTrigEffClosure.cxx`. Kept in lockstep by construction since 2026-08-18 (both read `ParamsSet`), but a signal-region change makes the MC-closure "data-like" variant STALE |
| `plotting_codes/single_b_analysis/plot_dr_vs_pair_pt_diagnostic.cxx` | :134, :137 | pp24 reco diagnostic -- still on `q*eta < 2.2` |
| `plot_reco_distr_singleb_vs_op_pp24.C` | :54-55 | pp24 reco diagnostic -- still on `q*eta < 2.2` |

> **Sites that carry ONLY the gap cuts** (single-muon windows + the pair-level `|eta^pair| < 2.2`
> added 2026-09-07), i.e. not the minv/pair-pT signal cuts, but which go stale with them and are
> easy to miss because they are not "the signal region": `RDFBasedHistFillingPP.cxx`
> `FillHistogramsGeneric` (generic + MC-data-comparison histograms) and its three template blocks;
> `RDFBasedHistFillingPythiaFullsim.cxx` `gap_truth` / `gap_reco` (which also build the
> `_gapcut_truth` generic family); `RDFBasedHistFillingPowhegTruth.cxx` `gap_truth`;
> `RDFBasedHistFillingPowhegFullsim.cxx` `gap_truth` (the `_single_b_pass_signal_truth_gapcut`
> family);
> `RDFBasedHistFilling/FillMCTrigEffHists.cxx` `kGapLeg` (Steps 2 & 4) and `kGapPair` (Step 3 +
> kn-split), the latter now sourced from `Utilities/MCTrigEffPairSelection.h`. Prefer grepping
> `FiducialGapCutExpr` and `PairFiducialEtaCutExpr` over trusting any line list here.

**Legacy / retired (leave, or fix for hygiene only):**
`SingleBAnalysis/SingleBAnalysisBase.cxx:25`,
`RDFBasedHistFilling/FillHistogramsCrossx_PP_clean.cxx:8` — both still use the
old `pair_eta < 2.2` form; not in the active chain.

> **⚠ ΔR-axis floors (only relevant when changing/removing the ΔR cut).** The
> ΔR-binned 2D histograms are booked with a **uniform axis from 0.05 to 1.0**:
> PP `:457/:519/:534`, PbPb `:1000/:1073/:1089/:1136/:1154`
> (`make_unif_edges(50, 0.05, 1.0)`). If `dr > 0.05` is removed, extend these
> axes to **0.0** or the newly-admitted small-ΔR pairs fall in underflow and the
> 2D ΔR distributions stay clipped. (Not needed for non-ΔR cut changes.)

---

## 3. [RERUN HIST FILLING] — order, scripts, outputs

**A. Data crossx (always required):**
1. `RDFBasedHistFilling/run_crossx_hist_filling_pp24.sh`
   → `dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root`
2. `run_crossx_hist_filling_pbpb23.sh`, `…pbpb24.sh`, `…pbpb25.sh`
   → `dimuon_data/pbpb_20YY/histograms_real_pairs_pbpb_20YY_single_mu4_no_trg_plots_nominal.root`
   (all three years — crossx & R_AA are always year-combined)

These regenerate the OS and SS signal-region histograms
(`h2d_crossx_…_w_signal_cuts` + `h2d_ss_…` for pp; `h3d_op/ss_…_vs_centr…` for
PbPb) that feed crossx plots **and** R_AA.

**B. Truth signal acceptance (required if acceptance/α is reported):**
3. `run_rdf_pythia_truth.sh` (modes: private, nonprivate 5.36, nonprivate 5.02)
   → `histograms_pythia_*_no_data_resonance_cuts.root` (+ `_near_away_divided`)
   — regenerates `h2d_sig_accept_{num,denom}_pt_eta` (and `_pt_150_eta`).

**C. MC fullsim reco-eff / det-response — code-consistency now, numerically
inert until real MC (see §4):**
4. `ENABLE_MC_TRIG_EFF=1 pipeline_pythia_fullsim_pp.sh full` — **pp24 reco-eff and detector
   response. RESULT-AFFECTING for pp, not optional** (§4 says so in as many words and this list
   used to omit it): the pp fullsim pair eps_reco has been the nominal correction since
   2026-08-18. Local pass over the LGD symlink farm, ~6-8 h, NOT Condor.
5. `plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C+(true)` →
   `pair_reco_eff_pp24_full.root`. **Belongs to no pipeline** — it must be run by hand after
   step 4, or `Utilities/PairRecoEffEvaluator.h` keeps applying the previous region's efficiency.
6. `pipeline_pythia_fullsim_overlay.sh` (hijing / zmumu / data) — PbPb reco-eff.
7. `pipeline_powheg_fullsim_single_muon.sh` — **single-muon, NO rerun** (does
   not use the pair signal cut).

> **Boundary warning (added 2026-09-08).** §1's "the NTuple stage is unchanged" premise does NOT
> hold for a change to the **muon pT threshold**, which is an NTuple-processing cut
> (`DimuonDataAlgCoreT.c` for data; the truth-pT gates in `PythiaAlgCoreT.c` / `PowhegAlgCoreT.c`
> and the reco/truth gates in `PythiaFullSimExtras.c` / `PowhegFullSimExtras.c` for MC). Such a
> change forces a FULL NTuple rerun in data and MC, and it drags POWHEG fullsim back into the
> blast radius — the carve-out in `pp24_all_vertex_pairs.md` rested on that change being
> reco-only, and a truth-pT cut is not. See `docs/tracking/mu_pt45_gap125_pairpt9_adoption.md`.

---

## 4. Reco-eff / det-response nuance (read before rerunning §3.C)

**pp24 (since 2026-08-18): the fullsim pair eps_reco IS the nominal.** It is
eps_reco(pair pT, pair eta, dR), built by
`plotting_codes/reco_effcy/build_pp24_fullsim_pair_reco_eff.C` on the canonical coarse axes and
read by `Utilities/PairRecoEffEvaluator.h`. For pp, §3.C is therefore RESULT-AFFECTING: a signal
cut change must rerun the fullsim reco-eff AND rebuild that file, or the correction and the data
describe different regions.

**Pb+Pb still applies the placeholder** (Run 2 TF1, pT & q*eta, dR-independent) -- **not** the
fullsim/overlay pair eps_reco. So, for Pb+Pb only:

- Removing/altering the signal cut changes the **set of pairs filled** in data
  crossx and the **truth acceptance** → §3.A and §3.B genuinely change results.
- The fullsim/overlay `pass_signal_*` definitions (§2) must be edited for
  consistency, but **re-deriving them does not change today's nominal crossx/
  R_AA** (placeholder is used). Treat §3.C as *bookkeeping until the real Run 3
  MC pair ε_reco(pair pT, pair η, ΔR) and det-response replace the placeholder*,
  at which point §3.C becomes result-affecting. Flagged in `docs/placeholder.md`.

---

## 5. [RERUN PLOTTING / FITTING] — consumers of the above

**Data crossx (pipeline stages, `pipelines/pipeline_{pp,pbpb}_crossx.sh`):**
- `plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx`
  → `plots/single_b_analysis/pp24/pp24_crossx_*.png`
- `plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx` (combined yrs)
  → `plots/single_b_analysis/pbpb_23_24_25_combined*/{TAA_weighted,counts}/*.png`
- `plotting_codes/single_b_analysis/plot_crossx_trig_corr_sanity.C` → sanity plots
- `plotting_codes/single_b_analysis/plot_crossx_reco_eff_stages.C`
  → `plots/sanity_check_crossx/` (raw/unfolded/reco/reco+trig stage overlays)

**R_AA (run manually, not in pipeline):**
- `RAA_plotting.cxx` (mode 6: PbPb 23+24+25+26 vs pp24, OS−SS) → `plots/single_b_analysis/RAA/*.png`

**MC–data comparison (PP crossx pipeline stage 8, optional):**
- `plotting_codes/mc_data_compr/plot_mc_data_compr.cxx`
  (+ `plot_mc_data_2D_hists_and_1D_proj.cxx`) → `plots/mc_data_compr/*.png`

**Signal acceptance / cutflow:**
- `plot_sig_accept_cutflow_above_60GeV.cxx` → `plots/single_b_analysis/{pythia,powheg}/*_sig_accept_above_60GeV*.png`
- `plot_signal_acceptance_pythia.cxx`, `plot_signal_acceptance_powheg.cxx`
  (`SignalAcceptancePlotter`) → 2D α(pair pT, pair η) + α-vs-pT-by-η subplots
  in `plots/single_b_analysis/{pythia,powheg}/`
  *(the 2D α(pair pT, pair η) plot itself is on the near-term to-add list.)*

**Reco-eff / det-resp plots (only meaningful once §3.C is result-affecting):**
- `plotting_codes/reco_effcy/PythiaFullsimRecoEffPlotter.cxx`,
  `plot_reco_effcy_pythia_fullsim_pp24.cxx`, det-response plotters.

**Internal note:** any synced figure embedding the above (crossx, R_AA,
acceptance) goes stale → re-run `/sync-note-figures` then `/check-note-sync`
(Gate G4). Figures live in `IntNotes/figures/` via the provenance manifest.

---

## 6. Ordered execution checklist (copy per change)

```
[ ] Edit cut in all §2 files (+ cutflow kCuts; + ΔR-axis floors if ΔR changes)
[ ] /review-analysis-code on the edits (quote analysis_overview §2 signal region)
[ ] Recompile (ACLiC) PP, PbPb, PythiaTruth, fullsim/overlay classes
[ ] *** IS THIS AN NTUPLE-STAGE CUT? *** (muon pT is; pair pT / gap window / minv are not)
[ ]   if yes: rerun NTuple processing FIRST -- data pp24 (both modes) + pbpb 23/24/25,
[ ]   Pythia truth, Pythia fullsim pp, HIJING overlay, POWHEG truth AND POWHEG fullsim.
[ ]   Nothing below is valid until it completes; see the §3.C boundary warning.
[ ] Rerun data crossx: pp24 + pbpb23/24/25            (§3.A)
[ ] Rerun pythia truth acceptance (3 modes)           (§3.B)
[ ] Rerun pp Pythia fullsim reco-eff + det-response    (§3.C/§4) -- RESULT-AFFECTING for pp
[ ] Rebuild pair_reco_eff_pp24_full.root by hand       (§3.C/§4) -- in no pipeline
[ ] (Real MC only) rerun fullsim/overlay reco-eff      (§3.C/§4)
[ ] Replot: pp crossx, pbpb crossx, sanity, stages    (§5)
[ ] Rerun RAA_plotting (mode 6)                        (§5)
[ ] (optional) MC–data comparison                      (§5)
[ ] Replot signal acceptance + cutflow                 (§5)
[ ] /review-plot on regenerated plots
[ ] /sync-note-figures + /check-note-sync              (§5)
[ ] Update docs (analysis_overview §2 if region changed; placeholder.md; roadmap)
```

## 7. Use for systematics

For a selection **systematic**, run §2–§6 with the varied cut into a **separate
output tag / backup dir** (do not overwrite nominal), then compare the varied
crossx/R_AA against nominal. The same blast radius applies; only data crossx
(§3.A) + truth acceptance (§3.B) + their plots are needed for a pure
acceptance/selection systematic (reco-eff placeholder is ΔR-independent, so it
cancels in the variation today).
