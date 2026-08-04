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
| 3 | Regenerate the single-muon trees with the current nominal selection (one-sided Δp/p) | IN PROGRESS |
| 4 | Plot the q·η spectra (pp24; PbPb combined ctr-integrated; PbPb combined ctr-binned) | TODO |
| 5 | `/review-plot` | TODO |

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

| Quantity | Value | Where |
|---|---|---|
| `eta_gap_cut1` | `0.135` — reject `fabs(eta) < 0.135` (charge-independent) | `ParamsSet.h:127` |
| `charge_eta_gap_cuts` | `{{0.56, 0.67}, {1.064, 1.29}, {-1.29, -1.12}}` — reject `charge*eta` inside any interval | `ParamsSet.h:145`, defined `:320` |
| `eta_gap_cut2` (COMMENTED OUT) | `{1.05, 1.29}`, applied only for `pt < 6` | `ParamsSet.h:128`; call site commented at `RDFBasedHistFillingData.cxx:149` |

Applied by `PassSingleMuonGapCut(eta, pt, charge)` →
`MuPairPassGapCut` = both muons pass (`RDFBasedHistFillingData.cxx:141–154`). Note the
`charge_eta_gap_cuts` loop is inside an `if` block — verify the guard before reuse.

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
