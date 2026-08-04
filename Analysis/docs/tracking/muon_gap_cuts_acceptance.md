# Muon detector-gap fiducial cuts + acceptance efficiency

Mode: **INVESTIGATION / decision support** (2026-08-04). Will be superseded by an
*implementation* doc (with a Physics Procedure) once the user finalises the gap-cut set.

## Objective

Provide the two inputs the user needs to finalise a **detector-gap fiducial cut** on single
muons, which would then be compensated by an **acceptance efficiency** ε_acc inserted into the
efficiency-correction chain:

1. The **single-muon q·η spectrum** after the nominal NTuple-processing selection, for
   - pp24 data,
   - PbPb 2023+2024+2025 data **combined**, centrality-integrated,
   - PbPb combined, centrality-binned (0–10 combined, 10–20, 20–30, 30–50, 50–80 %).
2. A complete **inventory of the gap cuts that already exist in the code base** (current RDF
   code + legacy pre-RDF code), with exact values and where each is (or is not) applied.

The user finalises the cut set from (1) + the data-driven mu4 trigger efficiencies. This doc
does NOT choose the cuts.

## Context

- Gap regions in the ATLAS muon spectrometer: the barrel crack at η≈0 (services/feet at
  z=0), the barrel↔endcap transition ≈1.05–1.3, and support-structure ("feet") regions. In
  these regions ε_reco and ε_trig are small and vary fast, so an efficiency correction there
  carries a large systematic; a fiducial cut + acceptance factor is the standard alternative.
- The analysis already has an unused/partially-used gap-cut machinery (see Findings F2).
- Single-muon q·η is already the canonical axis of the data-driven mu4 trigger efficiency, so
  a gap cut expressed in q·η composes directly with the existing efficiency chain.

## Sub-steps

| # | Step | Status |
|---|---|---|
| 1 | Locate/validate the single-muon trees; check staleness vs current nominal selection | DONE (F1, F3) |
| 2 | Inventory existing gap cuts | DONE (F2) |
| 3 | Regenerate the single-muon trees with the current nominal selection (one-sided Δp/p) | DONE (F5) |
| 4 | Plot the q·η spectra (pp24; PbPb combined ctr-integrated; PbPb combined ctr-binned) | DONE (F6) |
| 5 | `/review-plot` | IN PROGRESS |

## Accumulated Findings (append-only)

### F1 — Single-muon trees exist for all four samples (2026-08-04)

Produced 2026-06-23 by the four `NTupleProcessingCode/run_*_output_single_muon_tree.{sh,sub}`
scripts (commits `93761e0`, `47f81f0`):

| Sample | Path (dir `~/usatlasdata/dimuon_data/<dir>/`) | parts |
|---|---|---|
| pp24 | `pp_2024/single_muon_trees_pp_2024_part<N>_2mu4_mindR_0_02.root` | 12 |
| PbPb23 | `pbpb_2023/single_muon_trees_pbpb_2023_part<N>_single_mu4_mindR_0_02.root` | 4 |
| PbPb24 | `pbpb_2024/single_muon_trees_pbpb_2024_part<N>_single_mu4_mindR_0_02.root` | 2 |
| PbPb25 | `pbpb_2025/single_muon_trees_pbpb_2025_part<N>_single_mu4_mindR_0_02.root` | 6 |

Tree `muon_tree` (+ `muon_tree_ctr<1..20>` 5%-wide centrality trees for PbPb, written because
`turn_on_ctr_binned_tree_writing = true`). Branches: `ind, ev_num, pt, eta, phi, charge,
dP_overP, z0, d0, quality, passmu4, passmu4noL1, pass_tight, trk_*, ev_centrality, ev_FCal_Et`.
The `*_hadd.root` files are 599-byte failed hadds — **use a TChain over the parts**.

Provenance (`DimuonDataAlgCoreT.c:805`): a muon enters the tree if it is a member of a pair
that survives the **full pair-level nominal chain** — trigger match, both muons Combined+
Medium+IDCuts+MuonCuts, |η|<2.4, pT>4 GeV, Δp/p, |d0|/|z0 sinθ|, photoproduction veto (PbPb),
resonance veto (`resonance_cut_mode=1`) — de-duplicated by muon index. So it is *not* an
unbiased single-muon spectrum: it is the single-muon projection of the selected dimuon sample,
and it is trigger-biased (pp: 2mu4; PbPb: mu4). That is the right object here, because the gap
cut will be applied to exactly these muons.

`pass_tight` is never set on the data path (only in the fullsim `*FullSimExtras.c`), so the
working point must be applied at plot level via `quality & 16` (Tight) / `quality & 8`
(Medium). Trees are written at **Medium** (`requireTight = false` default).

### F2 — Existing gap cuts in the code base (2026-08-04)

**Single source of truth: `MuonObjectsParamsAndHelpers/ParamsSet.h`.**

| Quantity | Value | Applies to | Where |
|---|---|---|---|
| `eta_gap_cut1` | `0.135` — reject `fabs(eta) < 0.135` (charge-independent) | **all pT** | `ParamsSet.h:127` |
| `charge_eta_gap_cuts` | `{{0.56, 0.67}, {1.064, 1.29}, {-1.29, -1.12}}` — reject `charge*eta` inside any interval | **only `pt < 6` GeV** | `ParamsSet.h:145`, defined `:320` |
| `eta_gap_cut2` (COMMENTED OUT) | `{1.05, 1.29}` on `fabs(eta)`, also `pt < 6` | — | `ParamsSet.h:128`; call site commented at `RDFBasedHistFillingData.cxx:149` |

Applied by `PassSingleMuonGapCut(eta, pt, charge)` →
`MuPairPassGapCut` = both muons pass (`RDFBasedHistFillingData.cxx:141–154`).
**The pT gate is easy to miss:** the `charge_eta_gap_cuts` loop sits inside
`if (mpt < 6) { ... }` (`:144`), so above 6 GeV the ONLY gap rejection is `|η| < 0.135`.
Any reuse of these windows as a *fiducial* cut has to decide explicitly whether the pT
gate is kept — a pT-dependent fiducial region does not factorise into an acceptance
that is a pure function of q·η.

Note `|q·η| ≡ |η|`, so `eta_gap_cut1` is exactly the q·η window (−0.135, +0.135) and all
four existing rejections can be drawn on a single q·η axis.

Usage today: the gap cut is **NOT part of the signal selection**. It only produces an extra
diagnostic histogram family with the `_wgapcut` suffix (`RDFBasedHistFillingData.cxx:188–209`,
vars `Dphi`, `DR`, `DR_zoomin`), against the `_nogapcut` nominal. Labels
`ParamsSet::gapcut_labels = {"_nogapcut", "_wgapcut"}`.

Legacy pre-RDF code uses **the same `ParamsSet` values**, not different ones —
`OldMuonPairHistFillingCodePreRDF/{MuonPairPlottingPP.c:755, MuonPairPlottingPythia.c:341,
original/MuonPairPlotting.c:34, original/MuonPairPlottingPowheg.c:43}`,
`OldRDFBasedMuonPairPlotting/*.c:10`, `old/MuonPairPlottingOLDPP.c:31`. The only genuinely
different legacy variant is `original/MuonPairPlotting_coarse_ctr_bins.c:59`, which used a
single scalar `pms.eta_gap_cut` (`fabs(eta) > eta_gap_cut` on both muons) — that member no
longer exists.

**Second, independent set of "gap" regions** — the holes in the trigger-efficiency fit binning
`CommonEffcyConfig::q_eta_proj_ranges_fine_excl_gap` (Run 3, `CommonEffcyConfig.h:15`):
excluded q·η intervals are **(−1.3, −0.9)**, **(−0.1, 0.1)**, **(1.0, 1.3)**. Muons there are
not thrown away; they fall back to the unfitted 2D trigger-efficiency ratio
(`RDFBasedHistFillingData.cxx:628,646`; `FillMCTrigEffHists.cxx:330,353`). These do NOT match
`charge_eta_gap_cuts` — three coexisting definitions of "the gap":

| region | `charge_eta_gap_cuts` | `..._fine_excl_gap` |
|---|---|---|
| crack at 0 | `fabs(eta) < 0.135` (η, not q·η) | q·η ∈ (−0.1, 0.1) |
| feet-like | q·η ∈ (0.56, 0.67) | — (not excluded) |
| barrel/endcap + | q·η ∈ (1.064, 1.29) | q·η ∈ (1.0, 1.3) |
| barrel/endcap − | q·η ∈ (−1.29, −1.12) | q·η ∈ (−1.3, −0.9) |

Reconciling these into ONE definition is part of the decision the user is about to make.

### F2b — This decision is coupled to the LIVE bug in `pp_trig_eff_highpt_jump.md`

That doc was **REOPENED 2026-08-04** and is blocked on a user decision whose option **(d)**
is *"exclude the q·η gap from the fiducial region entirely"* — i.e. exactly the fiducial-cut
approach being considered here. The live defect: `EvaluateSingleMuonEffcyPtFitted`
(`RDFBasedHistFillingData.cxx:625–655`) returns the `-1.0f` sentinel when a muon's q·η is in
a `q_eta_proj_ranges_fine_excl_gap` hole **and** the 2D fallback bin is empty (which happens
at high pT); that propagates to `w_trig = 0`, so the pair is **silently dropped from the
nominal crossx histograms** — corrected/uncorrected falls to 0.37 in pp24 (η^pair ≈ −1.2,
pair pT 84–100 GeV) and to exactly 0 in PbPb. A gap fiducial cut chosen here would resolve
that bug as a by-product, provided the cut covers the trig-eff fit holes
(q·η ∈ (−1.3,−0.9), (−0.1,0.1), (1.0,1.3)) and not only the narrower
`charge_eta_gap_cuts` windows. **The two decisions should be taken together.**

### F3 — The existing trees are STALE w.r.t. the current nominal Δp/p cut (2026-08-04)

Commit `71fcf1c` (2026-07-16) made Δp/p **one-sided** (`dP_overP < thrsh`; the negative tail is
kept — the historical `fabs()` was a bug). The trees date from 2026-06-23, i.e. they were built
with `fabs(dP_overP) < thrsh` and are missing the negative tail. This is exactly the
"data-side Δp/p blast radius still stale" item in `mc_trigger_efficiency.md` (Remaining Work
3b). It matters here rather than being a formality: Δp/p is the ID/MS momentum imbalance, whose
negative tail is populated by muons losing more energy than expected — an η-dependent effect
that is largest precisely in the poorly-instrumented gap regions. ⇒ regenerate before plotting.

The only other post-2026-06-23 NTuple-processing change touching data, `0453074` (signed
`trk_charge`), is inert here: `turn_on_track_charge = false` by default
(`DimuonDataAlgCoreT.h:345`).

### F4 — Side observations (NOT part of this task; flagged for the user)

- **Cross-year FCal scaling is documented but not implemented.** `fcal_scale_pbpb_2024.root`
  exists in the data dir and `docs/placeholder.md` / `systematic_uncertainties.md` /
  `analysis_roadmap_2026_06.md` describe applying multiplicative 24/25→23 FCal scale factors,
  but `git log -S fcal_scale` shows the string has **only ever appeared in docs, never in
  code**. PbPb24/25 centrality currently comes from the raw `centrality` branch (PbPb24) or
  from the 2023 Glauber thresholds on unscaled FCal (PbPb25 pair-level recalculation).
- **Muon-level `ev_centrality` bypasses the PbPb25 pair-level recalculation.**
  `PbPbExtras.c:88` sets `m1/m2.ev_centrality = centrality` (raw branch) for all Run-3 years;
  the `MuonPairPbPb::PairValueCalcHook` recalculation (`MuonPairPbPb.h:53`) applies to
  `avg_centrality` only. Checked directly: the raw branch is NOT all-zero in the current PbPb25
  skim (mean 15.5, range −1…84, 4.6 % at 0), so the muon-level value is usable, but it is not
  the same object the pair-level analysis uses for PbPb25.

### F5 — Trees regenerated; the Δp/p staleness was real and measurable (2026-08-04)

All 24 Condor jobs (clusters 202–205) completed; 12 pp + 4/2/6 PbPb parts rewritten. Old
trees preserved under `<year>/single_mu_tree_bak_20260804_pre_dpop_fix/`. Direct check of
part 1 of each sample, NEW vs OLD:

| sample | N (new) | N (old) | muons with Δp/p < −0.12 (new) | (old) |
|---|---|---|---|---|
| pp24 | 1 014 411 | 978 010 | 18 541 (1.83 %) | **0** |
| PbPb23 | 305 723 | 288 644 | 8 788 (2.87 %) | **0** |
| PbPb24 | 1 024 201 | 966 170 | 29 937 (2.92 %) | **0** |
| PbPb25 | 978 218 | 920 550 | 29 729 (3.04 %) | **0** |

`Δp/p > +0.12` is 0 in both, as it must be. So the old trees were exactly the two-sided
`fabs` selection and the regeneration recovered a 1.8 % (pp) / 2.9–3.0 % (PbPb) population.

### F6 — q·η spectra produced (2026-08-04)

Macro `plotting_codes/single_b_analysis/plot_muon_q_eta_spectrum.cxx` (ACLiC), PNGs in
`~/usatlasdata/dimuon_data/plots/single_b_analysis/muon_gap_cuts/`:
`muon_q_eta_spectrum_pp24.png`, `..._pbpb_combined.png`,
`..._pbpb_combined_ctr_binned.png`, `..._pp_vs_pbpb_shape.png` (Medium WP → `medium/`).

Statistics at Tight WP: pp24 9 107 895 muons; PbPb 23+24+25 10 266 164 (2 767 001 +
1 953 036 + 5 546 127 before the 0–80 % centrality requirement).

Observed structure (both samples, all centralities):
- **A deep, narrow crack at q·η ≈ 0**, core roughly |q·η| ≲ 0.05, i.e. **substantially
  narrower than the existing `eta_gap_cut1 = 0.135` window**. Charge-symmetric.
- **A dip at q·η ≈ −1.15**, matching the existing `charge_eta_gap_cuts` window
  (−1.29, −1.12) well. It is **much deeper in pp (2mu4) than in PbPb (mu4)** and much
  deeper for pT < 6 GeV than above.
- **It appears at the same q·η for BOTH charges** — see F7 for the decisive quantitative
  test. It is therefore bending-direction dependent, not a fixed-|η| geometric hole.
- The positive-side window (1.064, 1.29) sits on much shallower structure than its
  negative-side counterpart.
- The (0.56, 0.67) "feet" window sits on weak structure in both samples.
- Sharp fall-off for |q·η| > 2.3 (edge of acceptance), strongest in pp.
- Centrality dependence of the SHAPE is weak (see the unit-area overlay panel): the
  fraction failing the existing gap cut runs 0.1022 → 0.0956 from 0–10 % to 50–80 %.

Fraction of single muons failing the EXISTING `PassSingleMuonGapCut` (counted per muon,
pT gate included): **pp24 7.39 %**, **PbPb 23+24+25 10.15 %** (0–10 % 10.22 %, 10–20 %
10.20 %, 20–30 % 10.11 %, 30–50 % 9.99 %, 50–80 % 9.56 %). The pp/PbPb difference comes
from the different trigger (2mu4 product vs mu4) and pT spectra, not from the cut.

### F7 — DECISIVE: the two features have different symmetry, so they need different cut variables (2026-08-04)

The "by muon charge" panel is suggestive but not conclusive. The decisive test histograms
**η** (not q·η) separately for μ⁺ and μ⁻ and measures the deficit depth
(`1 − N_core/N_sideband`, core |Δη|<0.05, side-bands 0.15–0.35 either side; uniform 0.02
η bins — a numeric diagnostic only, no plot, no new plotting binning):

| sample | charge | depth at η = −1.15 | depth at η = +1.15 | depth at η = 0 |
|---|---|---|---|---|
| pp24 | μ⁺ | **0.566** | 0.060 | 0.648 |
| pp24 | μ⁻ | 0.046 | **0.549** | 0.685 |
| PbPb 23+24+25 | μ⁺ | **0.283** | 0.084 | 0.536 |
| PbPb 23+24+25 | μ⁻ | 0.070 | **0.276** | 0.587 |

Interpretation — a **fixed-|η| geometric** gap must dip at BOTH η = ±1.15 for BOTH charges;
a **bending-direction** effect dips at η = −1.15 for μ⁺ and η = +1.15 for μ⁻ only (both are
q·η = −1.15). The measurement is unambiguous:

1. **The |η| ≈ 1.15 feature is a q·η effect, not a geometric gap.** Each charge dips on one
   side only, by 0.55 (pp) / 0.28 (PbPb), while the mirror side shows 0.05–0.08 — a factor
   ~7–9 asymmetry, far outside statistics (N ≈ 4.5–5.1 M per charge). ⇒ **a fiducial cut on
   this feature must be defined in q·η.** Cutting it in |η| would remove roughly twice the
   phase space needed, discarding the fully efficient mirror region. The existing
   `charge_eta_gap_cuts` are already correctly one-sided in q·η.
2. **The η ≈ 0 crack IS geometric and charge-symmetric** (0.648 vs 0.685 pp; 0.536 vs 0.587
   PbPb). ⇒ it is correctly cut in |η|, as `eta_gap_cut1` does.

So the two existing cuts use the right variable **each**, and the apparent inconsistency
noted in F2 (one in |η|, one in q·η) is physically justified, not an oversight.

Consistent with the KB: `physics/detector/atlas_run2_muon_trigger.md` records L1 barrel
coverage ≈80% with "gaps at η≈0, feet, support structures" — geometric and η-symmetric —
while the q·η dependence is the standard toroid bending-direction effect that motivates
ATLAS quoting muon trigger efficiency vs q·η in the first place. **RUN2-CROSSCHECK: no Run 2
figure of this raw single-muon q·η *yield spectrum* exists to compare against** (Run 2
reports *efficiencies* vs q·η, not the selected-pair muon spectrum), so no numerical Run 2
comparison is possible for the plotted quantity itself.

## Ruled Out (append-only)

- *Plotting from the existing 2026-06-23 trees* — stale on the one-sided Δp/p fix (F3), which
  is η-dependent and therefore not safely ignorable for a gap-region study.
- *Re-deriving the spectrum standalone from the raw NTUPs* — forbidden by the NTuple-Processing
  Provenance rule and unnecessary: the single-muon trees are exactly the intended output.

## Latest Stage

**Step 3 (in progress):** back up the 2026-06-23 single-muon trees, then resubmit all four
Condor productions (`run_{pp_24,pbpb_23,pbpb_24,pbpb_25}_output_single_muon_tree.sub`;
12/4/2/6 jobs) with the current nominal selection.

**Step 4 (next):** new macro `plotting_codes/single_b_analysis/plot_muon_q_eta_spectrum.C`.
Design decisions fixed up front:
- **q·η binning = `ParamsSet::makeEtaTrigEffcyBinning(1)`** — the existing canonical fine q·η
  axis (0.01 for |q·η|<0.2, 0.02 over [−1.3, 1.4] and [2.2, 2.4], 0.10 elsewhere), already used
  for the 2D trigger-efficiency histograms. NOT a new binning (§Binnings rule).
- Variable-width bins ⇒ plot dN/d(q·η) via `Scale(1., "width")`.
- **WP config var, default TIGHT** (`quality & 16`), Medium (`quality & 8`) reachable.
- Centrality bins from `RDFBasedHistFillingData::FindCtrSuffix` / `RDFBasedHistFillingPbPb.cxx:502`
  — `{0,5},{5,10},{10,20},{20,30},{30,50},{50,80}` with 0–5 and 5–10 merged per the request.
- PbPb 23+24+25 always combined (never per-year).
- Overlay the F2 gap regions as shaded bands so the user can judge them against the spectrum.
- No trigger-efficiency correction: the raw (trigger-biased) spectrum is what the user asked
  for, to be read together with the data-driven mu4 efficiencies.
