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
| 6 | TRUTH q*eta reference spectrum (MC, truth fiducial only) | DONE (F12) |
| 7 | Fix TLegend boxes covering curves in the new-gap-cut figures | DONE (F13) |
| 9 | Switch q·η spectra (data + truth) to linear y-scale; delete `old_gap_cuts/` (superseded by F10's settled cut) | DONE (F15) |

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

**Second, independent set of "gap" regions** — the q·η intervals with NO fitted turn-on, i.e.
the complement of `CommonEffcyConfig::q_eta_proj_ranges_fine_excl_gap` (Run 3,
`CommonEffcyConfig.h:15`): **(−1.3, −0.9)**, **(−0.1, 0.1)**, **(1.0, 1.3)** — and, easily
missed, **(2.2, 2.4)**, because the fitted set is *one-sided* (−2.4 ≤ q·η < 2.2) and simply
stops at 2.2. Muons in all four are not thrown away; they fall back to the unfitted 2D
trigger-efficiency ratio (`RDFBasedHistFillingData.cxx:628,646`;
`FillMCTrigEffHists.cxx:330,353`) — and, when that bin is empty, to the `-1.0f` sentinel of
F2b. These do NOT match `charge_eta_gap_cuts` — three coexisting definitions of "the gap":

| region | `charge_eta_gap_cuts` | no fitted trigger efficiency |
|---|---|---|
| crack at 0 | `fabs(eta) < 0.135` (η, not q·η) | q·η ∈ (−0.1, 0.1) |
| feet-like | q·η ∈ (0.56, 0.67) | — (covered by a fit) |
| barrel/endcap + | q·η ∈ (1.064, 1.29) | q·η ∈ (1.0, 1.3) |
| barrel/endcap − | q·η ∈ (−1.29, −1.12) | q·η ∈ (−1.3, −0.9) |
| forward + edge | — (not rejected) | q·η ∈ (2.2, 2.4) |

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
- **Muon-level `ev_centrality` bypasses the PbPb25 pair-level recalculation — but the bypass
  is a measured NO-OP.** `PbPbExtras.c:88` sets `m1/m2.ev_centrality = centrality` (raw
  branch) for all Run-3 years, whereas `MuonPairPbPb::PairValueCalcHook` (`MuonPairPbPb.h:53`)
  recalculates `avg_centrality` from FCal for year 25. Recomputing
  `GetCentralityPbPb2023(ev_FCal_Et)` per muon and comparing with the raw branch over the full
  Tight samples: **PbPb25 0 differences / 5 546 127; PbPb24 0 / 1 953 036; PbPb23 2 / 2 767 001**
  (max \|Δ\| = 1, i.e. float boundary ties). So the centrality axis of the plots here IS on the
  2023 Glauber thresholds for all three years, and the earlier concern in this bullet was
  overstated. (The comment at `MuonPairPbPb.h:51` that the PbPb25 `centrality` branch "is all
  zeros in the skim" is also stale for the current skim: mean 15.5, range −1…84, 4.6 % at 0.)

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
`..._pbpb_combined_ctr_binned.png`, `..._pp_vs_pbpb_shape.png`. The Medium-WP set is
generated by `plot_muon_q_eta_spectrum.cxx+(false)` into `muon_gap_cuts/medium/`.

**Working-point caveat.** The WP is applied PER MUON (`quality & 16`), which is not the
nominal pair-level requirement (`m1.quality & m2.quality & 16`,
`DimuonDataAlgCoreT.c:584`) — the single-muon tree is de-duplicated and carries no partner
information, so per-muon is the only form available, and ~5 % of these muons come from
pairs the nominal Tight analysis drops. That is the right population for choosing a
*single-muon* fiducial cut, but it is NOT the nominal pair selection. `use_tight_wp=false`
is effectively a **no-op** (verified: pp24 9 608 214 in tree → 9 608 214 plotted): the
trees are already written with pair-level Medium applied.

**Binning-resolution caveat.** Over |q·η| ∈ [1.4, 2.2] the canonical
`makeEtaTrigEffcyBinning` step is 0.10, coarser than the real structure there (in 0.02
bins pp runs 33 209 at +1.45 vs 48 160 at +1.63, a ~45 % modulation the 0.10 bins average
over). A cut boundary in that region cannot be read off these figures to better than 0.1.
The canonical binning is deliberately NOT changed to suit a plot.

Statistics at Tight WP: pp24 9 107 895 muons. PbPb 23+24+25: 10 266 164 before the 0–80 %
centrality requirement (2 767 001 + 1 953 036 + 5 546 127), **10 237 444 after** — the
latter is the plotted and counted population.

Observed structure (both samples, all centralities):
- **A deep, narrow crack at q·η ≈ 0**, core roughly |q·η| ≲ 0.05, i.e. **substantially
  narrower than the existing `eta_gap_cut1 = 0.135` window**. Its *depth* is
  charge-symmetric, but its *profile* is not — see F7 conclusion 2.
- **A dip near q·η ≈ −1.1**, **NOT well matched by the existing `charge_eta_gap_cuts`
  window (−1.29, −1.12)** — see the measured profile below. Much deeper in pp (2mu4) than
  in PbPb (mu4), and much deeper for pT < 6 GeV than above (visible directly in the
  pT-split ratio panel).
- **It appears at the same q·η for BOTH charges** — see F7 for the decisive quantitative
  test. It is therefore bending-direction dependent, not a fixed-|η| geometric hole.
- The positive-side window (1.064, 1.29) sits on much shallower structure than its
  negative-side counterpart.
- The (0.56, 0.67) "feet" window sits on weak structure in both samples.
- **The acceptance edge is ONE-SIDED in q·η, not symmetric in |q·η|** (corrected — the
  first version of this note wrongly wrote "|q·η| > 2.3"). Measured in 0.02-wide slices:

  (Tight WP; PbPb over the same gated 0–80 % population as everything else.)

  | slice | pp24 | PbPb 23+24+25, 0–80 % |
  |---|---|---|
  | N(−2.40, −2.38) | 18 568 | 28 443 |
  | N(−2.22, −2.20) | 25 998 | 39 380 |
  | N(+2.20, +2.22) | 19 808 | 29 520 |
  | N(+2.38, +2.40) | **1 737** | **15 412** |

  The negative side retains 71 % (pp) / 72 % (PbPb) of its 2.2 value out at 2.4; the
  positive side retains only **8.8 % (pp)** / 52 % (PbPb) — a factor **11.4** collapse in
  pp. Same bending-direction physics as the −1.1 dip. **A symmetric `|q·η| < 2.3` cut would
  discard a fully populated region near q·η ≈ −2.4.**
- Centrality dependence of the SHAPE is weak — the shape/0–10 % ratio panel stays within
  roughly ±10 % across the whole axis; the fraction failing the existing gap cut runs
  0.1022 → 0.0956 from 0–10 % to 50–80 %.

**Measured profile of the barrel/endcap dip vs the existing window.** pp24 Tight, 0.02-wide
q·η bins **labelled by LOW EDGE**, μ⁺+μ⁻; local plateau ≈ 48 527 per bin (mean of the
(−1.60,−1.50) and (−0.90,−0.80) side-bands). Second column = yield/plateau:

```
 [-1.30,-1.28): 45156  0.93        [-1.14,-1.12): 16836  0.35
 [-1.28,-1.26): 42945  0.88        [-1.12,-1.10): 15419  0.32  <-- MINIMUM
 [-1.26,-1.24): 39629  0.82        [-1.10,-1.08): 18143  0.37
 [-1.24,-1.22): 35846  0.74        [-1.08,-1.06): 25066  0.52
 [-1.22,-1.20): 32326  0.67        [-1.06,-1.04): 31141  0.64
 [-1.20,-1.18): 29067  0.60        [-1.04,-1.02): 36415  0.75
 [-1.18,-1.16): 24176  0.50        [-1.02,-1.00): 39910  0.82
 [-1.16,-1.14): 20029  0.41
```

The minimum sits in **[−1.12, −1.10)** at 32 % of plateau; the half-depth (50 %) points are
≈ −1.17 and ≈ −1.08, so the **measured half-depth window is ≈ (−1.17, −1.08), centre
≈ −1.125**. The existing window (−1.29, −1.12) has centre −1.205 — **offset by ≈0.08**, and
mismatched at both ends:
- it **includes** −1.29…−1.22, where the yield is only **7–33 % depleted** (0.93 → 0.67 of
  plateau) — i.e. it cuts away largely good acceptance;
- it **excludes** −1.12…−1.00, which is still **68 % → 18 % depleted** (0.32 → 0.82) —
  including the minimum bin itself, since the window's upper edge −1.12 is exactly where
  the minimum bin begins.

In PbPb the same dip is far shallower and slightly displaced: minimum 35 603 in
[−1.18, −1.16) against a plateau of 53 187, i.e. **67 % of plateau rather than 32 %**.

Fraction of single muons failing the EXISTING `PassSingleMuonGapCut` (counted per muon over
exactly the plotted population, pT gate included): **pp24 7.39 %** (672 686 / 9 107 895),
**PbPb 23+24+25 0–80 % 10.16 %** (1 039 883 / 10 237 444); by centrality 0–10 % 10.22 %,
10–20 % 10.20 %, 20–30 % 10.11 %, 30–50 % 9.99 %, 50–80 % 9.56 %. The pp/PbPb difference
comes from the different trigger (2mu4 product vs mu4) and pT spectra, not from the cut.

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
2. **The η ≈ 0 crack affects BOTH charges — but its profile is NOT |η|-symmetric, so the
   right cut variable for it is NOT yet settled.** (Corrected: the first version of this
   note concluded "geometric and charge-symmetric ⇒ correctly cut in |η|". That conclusion
   rested on the *core depth* statistic above, which is symmetric by construction — it
   integrates |Δη| < 0.05 — and is therefore blind to skew.) Both charges do dip at η ≈ 0
   (0.648 vs 0.685 pp; 0.536 vs 0.587 PbPb), so a geometric component is real. But the
   q·η profile of the crack is strongly skewed, in the SAME direction for both charges
   (0.005-wide q·η bins, Tight WP):

   | | pp24 | PbPb 0–80 % |
   |---|---|---|
   | N(−0.050, −0.045) | 3 040 | 4 189 |
   | N(+0.045, +0.050) | 7 415 | 8 442 |
   | ratio (+/−) | **2.44** | **2.02** |
   | minimum at q·η | −0.0225 | −0.0275 |
   | plateau (±0.19) | ≈ 13 600 | ≈ 13 100 |

   Split by charge, in the SAME 0.005-wide bins (these sum exactly to the combined row —
   an earlier version of this table quoted 0.01-wide values here and was ~2× too large):

   | | N(−0.050,−0.045) | N(+0.045,+0.050) | ratio (+/−) |
   |---|---|---|---|
   | pp24 μ⁺ | 1 512 | 4 041 | 2.67 |
   | pp24 μ⁻ | 1 528 | 3 374 | 2.21 |
   | PbPb μ⁺ | 2 147 | 4 572 | 2.13 |
   | PbPb μ⁻ | 2 042 | 3 870 | 1.90 |

   **Both** charges are depleted on the q·η < 0 side — a bending-direction signature, not
   an η-symmetric one (the two charges do differ somewhat, 2.67 vs 2.21 in pp, so the
   effect is not purely a function of q·η either). So the crack carries a q·η component
   superimposed on a geometric gap, its minimum is offset from 0 by ≈0.02, and a symmetric
   `|η| < 0.135` window is not obviously the right shape. **This is a live question for the
   user, not a settled point.**

So the |η| ≈ 1.15 feature definitely needs a q·η cut (conclusion 1). Whether the crack is
best cut in |η| or in q·η is open; what F2 flagged as an apparent inconsistency between the
two existing cut variables is at least *partly* justified, but not fully established.

Consistent with the KB: `physics/detector/atlas_run2_muon_trigger.md` records L1 barrel
coverage ≈80% with "gaps at η≈0, feet, support structures" — geometric and η-symmetric —
while the q·η dependence is the standard toroid bending-direction effect that motivates
ATLAS quoting muon trigger efficiency vs q·η in the first place. **RUN2-CROSSCHECK: no Run 2
figure of this raw single-muon q·η *yield spectrum* exists to compare against** (Run 2
reports *efficiencies* vs q·η, not the selected-pair muon spectrum), so no numerical Run 2
comparison is possible for the plotted quantity itself.

### F8 — PROVISIONAL cut set written to ParamsSet.h (2026-08-04)

At the user's request, a proposed default is now declared in
`MuonObjectsParamsAndHelpers/ParamsSet.h` (declaration ~line 146, definition ~line 358):

```cpp
static std::vector<std::pair<float,float>> single_mu_fiducial_gap_cuts;
// = {{-1.20f, -1.05f}, {-0.06f, 0.06f}, {2.20f, 2.40f}}
```

**Nothing reads it yet — adding it changes no result.** `charge_eta_gap_cuts` and
`eta_gap_cut1` are untouched and `PassSingleMuonGapCut` still behaves exactly as before.
The user will fine-tune before it is wired into any selection.

Rationale per window is in the header comment; the measurements are F6/F7. Measured cost at
Tight WP:

| forward window | pp muons | pp pairs | PbPb muons | PbPb pairs |
|---|---|---|---|---|
| `{2.20, 2.40}` (as written) | 4.28 % | 8.38 % | 6.45 % | 12.49 % |
| `{2.30, 2.40}` (looser alternative) | 3.29 % | 6.47 % | 5.13 % | 10.00 % |

against 7.39 % / 10.16 % of muons for the current `PassSingleMuonGapCut`. The forward slice
is **already outside the signal region today**, so adopting `{2.20, 2.40}` costs nothing
relative to the current selection; loosening to `{2.30, 2.40}` would *recover* 1.36 % (pp) /
2.24 % (PbPb) of muons. That is the choice the user flagged as still open.

**Intended scope (user, 2026-08-04; NOT implemented).** This vector is to become the single
source of the gap/fiducial definition — applied in the signal selection **and** in the trigger-
and reconstruction-efficiency evaluations — **replacing** the standalone per-muon `q·η < 2.2`
cut currently repeated across `RDFBasedHistFilling{PP,PbPb,PythiaTruth,PythiaFullsim,
PythiaFullsimOverlay,PowhegTruth,PowhegFullsim}.cxx` (13 call sites listed below). Because the
`{2.20, 2.40}` window reproduces that cut exactly (muons already satisfy |η| < 2.4), the
replacement is **yield-neutral by construction** at the current edge — which makes it a clean
swap to validate before any edge is moved. Wiring it in is a signal-selection change with the
full rerun blast radius in `signal_selection_change_impact.md`.

**Correction to the earlier F6 forward-edge discussion.** A per-muon signal cut
`m1.charge*m1.eta < 2.2 && m2.charge*m2.eta < 2.2` **already exists** in every pipeline —
`RDFBasedHistFillingPP.cxx:421,514,684`, `RDFBasedHistFillingPbPb.cxx:964,1033,1323`,
`RDFBasedHistFillingPythiaTruth.cxx:471`, `RDFBasedHistFillingPythiaFullsim.cxx:137,139`,
`...FullsimOverlay.cxx:80,82`, `...PowhegTruth.cxx:158`, `...PowhegFullsim.cxx:185,186`.
Consequences:
1. A forward-edge fiducial window is **redundant against the current selection** — the
   one-sided collapse beyond q·η ≈ +2.3 is already outside the signal region. It is
   nevertheless **included** in the proposed vector as `{2.20, 2.40}` (user, 2026-08-04),
   because the vector is intended to *replace* that standalone cut rather than sit alongside
   it; writing it at 2.20 makes the swap yield-neutral, after which the edge can be moved.
2. The `(2.2, 2.4)` no-fitted-efficiency hole (F2) is likewise already outside the signal
   region, so it does **not** contribute to the F2b `w_trig = 0` bug in the cross-section.
   Only `(−1.3,−0.9)`, `(−0.1,0.1)` and `(1.0,1.3)` do.
3. The q·η > 2.2 muons visible in the plots are there because the single-muon trees carry the
   *pair-selection* cuts, not the crossx *signal* cuts (minv, pair pT, q·η) — the plotted
   population is deliberately upstream of the signal region.

**Recommendation recorded for the user's review (NOT applied):** the existing windows are not
adequate — (i) the (−1.29,−1.12) window is offset by ≈0.08, cutting acceptance that is only
7–33 % depleted while keeping the minimum bin; (ii) the pT<6 gate makes the fiducial region
pT-dependent so ε_acc does not factorise in q·η; (iii) they leave 12.7 % (pp) / 12.8 % (PbPb)
of muons with no fitted trigger efficiency, so they do not resolve F2b. **But cutting all the
no-fitted-efficiency regions is NOT the answer either**: it costs 17.5 % (pp) / 20.1 % (PbPb)
of muons ≈ **32 % / 36 % of pairs**. Those holes are a *fitting-configuration* choice, not a
measurement limitation (the 2D map covers them and only runs out of statistics at high pT), so
the right fix for F2b is to extend the fit binning — or option (a) nearest-fitted-bin — not to
redefine the acceptance around it. Keep the two decisions separate.

### F9 — Plot sets split old vs new gap cuts (2026-08-04)

`plot_muon_q_eta_spectrum.cxx` now takes a second argument, `use_new_gap_cuts`, selecting
WHICH candidate definition is overlaid and which one the quoted "fraction removed" refers to.
**The spectra themselves are identical in both cases** — only the overlay and the fractions
change. Output split accordingly:

```
muon_gap_cuts/old_gap_cuts/[medium/]   legacy PassSingleMuonGapCut  (red/orange bands)
muon_gap_cuts/new_gap_cuts/[medium/]   proposed fiducial windows    (green bands)
```

Both sets regenerated from the same macro build (4 PNGs x 2 WPs x 2 cut sets = 16).
Run as `plot_muon_q_eta_spectrum.cxx+(use_tight_wp, use_new_gap_cuts)`.

The proposed windows are **read from `ParamsSet::single_mu_fiducial_gap_cuts`** — nothing is
retyped or overridden in the macro, so re-tuning means editing `ParamsSet.h` alone and
re-running.

**Forward edge is now 2.30, not 2.20.** The value in `ParamsSet.h` was changed to
`{2.30f, 2.40f}` outside this session (it arrived in sibling commit `c426437`, a pair-pT
binning commit that also swept up the then-uncommitted `single_mu_fiducial_gap_cuts`
addition). That matches the user's request to plot the 2.3 variant, so the macro simply reads
it. **The header comment in `ParamsSet.h` was corrected accordingly** — it previously argued
the 2.20 edge made a future swap yield-neutral, which is no longer true.

**Consequence of 2.30 that must not be lost:** the vector is meant to REPLACE the standalone
per-muon `q·η < 2.2` signal cut. At 2.30 that swap is **not** yield-neutral — it recovers
q·η ∈ (2.2, 2.3), i.e. +1.36 % of pp and +2.24 % of PbPb single muons. **That recovered slice
has no fitted trigger efficiency today** (the fitted ranges in
`CommonEffcyConfig::q_eta_proj_ranges_fine_excl_gap` stop at 2.2), so adopting 2.30 REQUIRES
extending the fit binning to cover (2.2, 2.3) — otherwise those muons hit the `-1.0f`
sentinel and are silently dropped (F2b / `pp_trig_eff_highpt_jump.md`). Reverting the edge to
2.20 would instead make the swap exactly yield-neutral and need no fit change.

Measured fraction of single muons removed, per cut set (per muon; pairs assume both must pass):

| cut set | pp muons | pp pairs | PbPb 0–80 % muons | PbPb pairs |
|---|---|---|---|---|
| legacy `PassSingleMuonGapCut` | 7.39 % | 14.2 % | 10.16 % | 19.3 % |
| proposed, Tight | **3.29 %** | 6.5 % | **5.13 %** | 10.0 % |
| proposed, Medium | 3.36 % | 6.6 % | 5.72 % | 11.1 % |

Per centrality (proposed, Tight): 0–10 % 5.09 %, 10–20 % 5.17 %, 20–30 % 5.15 %,
30–50 % 5.16 %, 50–80 % 5.27 % — essentially flat, unlike the legacy set which drifts
0.1022 → 0.0956 over the same range.

### F10 — USER DECISIONS: the cut set is SETTLED (2026-08-12)

The three open questions from F8/F9 were put to the user and answered. The fiducial definition
is now fixed and matches what is already in `ParamsSet.h` — **no code change was needed**:

```cpp
single_mu_fiducial_gap_cuts = {{-1.20f, -1.05f}, {-0.06f, 0.06f}, {2.30f, 2.40f}};
```

1. **Forward edge = 2.30** (not 2.20). The user chose the looser edge, accepting that the swap
   for the standalone per-muon `q·η < 2.2` signal cut is therefore **not** yield-neutral: it
   recovers q·η ∈ (2.2, 2.3), +1.36 % of pp and +2.24 % of PbPb single muons.
   **Consequence, now a hard PREREQUISITE for wiring in:** that recovered slice has no fitted
   trigger efficiency (`CommonEffcyConfig::q_eta_proj_ranges_fine_excl_gap` stops at 2.2), so
   the fit binning MUST be extended to cover (2.2, 2.3) *before* the vector replaces the
   standalone cut — otherwise those muons hit the `-1.0f` sentinel and are silently dropped
   (F2b / `pp_trig_eff_highpt_jump.md`). The user did NOT ask for that extension yet; it is
   deferred with the wiring-in, not waived.
2. **Crack window stays in q·η**, `{-0.06, +0.06}`, not switched to |η|. One variable for all
   three windows; the measured skew (minimum at q·η ≈ −0.02, ratio 2.44 pp / 2.02 PbPb at
   ∓0.05, F7) is accepted rather than modelled with a per-window variable flag.
3. **Wiring in is HELD.** The replacement at the 13 call sites + the trigger/reco-efficiency
   evaluations is not started and needs an explicit go-ahead. Nothing reads
   `single_mu_fiducial_gap_cuts`, so no result has moved.

### F11 — Crack window widened to (−0.10, +0.06); dashed no-fit markers removed (2026-08-12)

**(a) Central crack window is now ASYMMETRIC.** `single_mu_fiducial_gap_cuts` central window
changed `{-0.06, +0.06}` → `{-0.10, +0.06}` (user, 2026-08-12). The full set is now

```cpp
{{-1.20f, -1.05f}, {-0.10f, 0.06f}, {2.30f, 2.40f}};
```

The user's instruction as first written was `(-1., -0.06)`. Measured before applying anything:
that window removes **23.81 % (pp) / 24.51 % (PbPb)** of single muons on its own — 25.99 % /
28.30 % for the whole set, i.e. **45.2 % / 48.6 % of PAIRS** — and it would have cut a region
that is *enriched*, not depleted. Slice-by-slice N(negative)/N(mirror positive) in 0.1 steps:

| slice | pp | PbPb |
|---|---|---|
| (−0.10, 0.00) | **0.509** | **0.578** |
| (−0.20, −0.10) | 0.885 | 0.872 |
| (−0.30, −0.20) … (−1.00, −0.90) | 1.03 – 1.46 | 1.08 – 1.39 |
| (−1.10, −1.00) | 0.772 | 1.059 |
| (−1.20, −1.10) | 0.599 | 0.932 |

Only (−0.10, 0.00) is depleted near zero — and only on the negative side, in the same sense
for both charges (F7). It would also have left a surviving 1.0 % (pp) / 1.3 % (PbPb) sliver at
[−1.05, −1.00] between the two windows. Put to the user with these numbers; the user chose
`(-0.10, +0.06)`, which covers the whole measured depletion and nothing else.

Cost at Tight WP, measured twice independently (a scratch scan over the trees, then the macro's
own fill-loop counter — they agree exactly):

| set | pp muons | pp pairs | PbPb muons | PbPb pairs |
|---|---|---|---|---|
| legacy `PassSingleMuonGapCut` | 7.39 % | 14.2 % | 10.16 % | 19.3 % |
| fiducial, symmetric `{-0.06,+0.06}` (previous) | 3.29 % | 6.47 % | 5.13 % | 10.00 % |
| **fiducial, `{-0.10,+0.06}` (current)** | **3.79 %** | **7.44 %** | **5.60 %** | **10.89 %** |

Per-window: crack 1.62 % (pp) / 1.81 % (PbPb), dip 1.80 % / 2.88 %, forward 0.37 % / 0.91 %.
Medium WP: 3.90 % (pp) / 6.40 % (PbPb) of muons. Per centrality (Tight) 5.55 → 5.76 % over
0–10 → 50–80 %, still essentially flat.

**(b) THE VECTOR IS NO LONGER INERT — this edit changes results.** F8/F10 said nothing read
`single_mu_fiducial_gap_cuts`. That is now WRONG: a sibling session wired it into the
**trigger-efficiency chain** on 2026-08-04 (`mc_trigger_efficiency.md`), and the stale
"nothing reads this vector" claim has been corrected in the `ParamsSet.h` header comment.
Live consumers:
- `RDFBasedHistFilling/FillMCTrigEffHists.cxx:641–663` — MC trig-eff, Steps 1–4, single /
  leg / pair forms via `ParamsSet::FiducialGapCutExpr`.
- `RDFBasedHistFillingPP.cxx:172`, `RDFBasedHistFillingPbPb.cxx:342` — data tag-and-probe,
  applied to the **probe only** (user decision; the tag is left uncut).
- `RDFBasedHistFillingData.cxx:72–74` — reads the forward window's lower edge to check it
  against the coarse q·η binning.
- `CommonEffcyConfig.h:48,63` — the coarse binning's top edge (2.30) must track the forward
  window or the evaluator throws. **The central window is not referenced there, so widening
  the crack does NOT break that consistency check** (the coarse binning has an edge at 0.0,
  and (−0.10, +0.06) sits inside the (−0.5, 0.0) and (0.0, +0.5) bins).

⇒ **The MC trigger efficiency and the data tag-and-probe efficiencies are now STALE** and must
be rerun. NOT done here: that is the sibling session's active pipeline, and rerunning it
unilaterally would collide with its in-flight work. Flagged to the user instead.

**(c) Dashed "no fitted trigger efficiency" markers removed from all 16 PNGs.** The trigger-
efficiency code now fits in the **contiguous** `CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap`
(10 bins, −2.4 ≤ q·η < 2.3, no holes, top edge tracking the fiducial cut), so the regions those
lines marked no longer exist — every surviving muon has a fitted turn-on and there is no
`-1.0f` sentinel to fall into. Verified in `CommonEffcyConfig.h:65–76` before removing. Deleted
from `plot_muon_q_eta_spectrum.cxx`: the `FineExclGapHoles()` helper, the dashed-`TLine` overlay
in `DrawGapMarkers`, the `"kept, but no fitted trigger efficiency"` `TLegend` entry, the console
report of the holes, and the now-unused `CommonEffcyConfig.h` include + `pairToSuffix` link shim.
Legend boxes shrunk by one entry's height and `CanvasGapKey` dropped from 2/3 to 1/2 columns so
the remaining entries keep their spacing. All 16 PNGs regenerated and visually checked.

### F12 — TRUTH q·η reference: the three gap features are 100 % detector, and ε_acc = 0.913 (2026-08-13)

New plot set `muon_gap_cuts/truth_q_eta/` from
`plotting_codes/single_b_analysis/plot_muon_truth_q_eta_spectrum.cxx`:
`muon_truth_q_eta_spectrum_pp24_fullsim.png`, `..._pbpb_overlay.png`,
`..._pp_vs_pbpb_shape.png`. Truth muons in the truth fiducial (`truth_pt > 4`,
`|truth_eta| < 2.4`) with **no reco match, no working point, no trigger** — the spectrum the
data figures of F6 cannot provide. Cross-section-weighted by `ev_weight`; same canonical
`makeEtaTrigEffcyBinning(1)` axis as F6, so truth and data panels are cell-for-cell
comparable. Sample scope was the user's call: pp24 primary **plus** the overlay; a
truth-vs-data comparison panel was offered and declined.

**Input provenance.** `muon_tree` in NOMINAL mode (`PythiaFullSimExtras.c:395–412`), whose
gate is exactly the truth fiducial with no reco match — it is the single-muon reco-efficiency
denominator. The `store_mc_trigger` twin of the same tree is NOT interchangeable (its gate is
`reco_match && pt > 3 && |eta| < 2.6`, which would put the reconstruction efficiency straight
back in). Truth-seeded from the Pythia-only block, so HIJING muons never enter.

| sample | truth muons | Σ weights (nb) |
|---|---|---|
| pp24 Pythia fullsim, FULL `_pdf` production | 17 105 648 | 168.3 |
| HIJING overlay, PbPb23 test sample, centrality 0–13 % | 107 786 | 144.1 |

**The overlay input was STALE and was regenerated here** (same failure mode as F3, caught
before plotting). The on-disk copy was dated 2026-06-10 and so predated `8c917b4`
(2026-07-14 — the Pythia truth index guard *leaked the HIJING generator block* when the
Pythia Geant4 block is empty, i.e. a change to which truth muons exist at all) and `216f750`
(2026-07-14, `isTestSample` / isospin, which sets `ev_weight`). Rerun with
`run_pythia_fullsim_overlay_single_muon.sh` → the new
`..._hijing_overlay_pbpb23_..._single_muon.root` (the old file also carried the legacy
`hijing_overlay_pp24` label). Measured effect: 107 796 → 107 786 muons, i.e. the HIJING leak
was **10 muons** in this sample — small, but it could not be known to be small without
rerunning. The pp24 FULL file (2026-07-20) postdates both fixes and is current: the only later
NTP commits are `c0043bd` (inside `if (store_mc_trigger)` only) and `f990173` (adds `n_vtx`).

**★ RESULT 1 — every feature the fiducial cut targets is purely detector.** The dip finder
(bins below 60 % of the local 21-bin median) reports **nothing anywhere on the axis, in either
sample**: the truth spectrum is a smooth, charge-symmetric arc. Measured against the same
diagnostic used on data in F6/F7 (0.02-wide q·η slices, cross-section weighted):

| feature | data pp24 (F6/F7) | **truth pp24** |
|---|---|---|
| barrel/endcap dip, minimum bin [−1.12, −1.10) | **0.32** of plateau | **1.036** |
| whole dip window [−1.20, −1.05], mean | ~0.4 of plateau | **1.034** |
| crack window [−0.10, +0.06], mean / side-band | deep and skewed (min ≈ 0.35, ratio 2.44 at ∓0.05) | **1.007** |
| forward N(+2.38,+2.40) / N(−2.40,−2.38) | **0.094** (factor 11.4 collapse) | **1.005** |

So the −1.11 dip, the η ≈ 0 crack **and** the one-sided forward collapse have **no generator-level
counterpart whatsoever** — they are entirely trigger × reconstruction. This is what makes the
truth spectrum usable as the ε_acc denominator, and it independently confirms the F7 reading
that these are detector/bending-direction effects rather than kinematics.

**★ RESULT 2 — the truth-level acceptance cost, ε_acc.** Fraction of truth muons removed by
`single_mu_fiducial_gap_cuts` (cross-section weighted; raw fractions agree to <0.4 %):

| window | pp24 | PbPb overlay |
|---|---|---|
| [−1.20, −1.05] (barrel/endcap dip) | 3.36 % | 3.32 % |
| [−0.10, +0.06] (crack) | 4.16 % | 4.15 % |
| [+2.30, +2.40] (forward edge) | 1.16 % | 1.32 % |
| **all windows** | **8.67 %** | **8.79 %** |
| **ε_acc = 1 − removed** | **0.9133** | **0.9121** |

**Two things to carry forward:**
1. **ε_acc must be built from TRUTH, not from the data spectra of F11a.** The same cut removes
   **8.67 % of truth muons but only 3.79 % of reconstructed pp muons** (5.60 % PbPb) — a factor
   **2.3 (pp) / 1.6 (PbPb)**. The data fractions are small precisely *because* the detector has
   already removed most of that population; using them as ε_acc would under-correct by that
   factor. The F11a numbers are the cost in *reconstructed yield*, which is a different (and
   also useful) quantity — the two must not be interchanged.
2. **ε_acc is essentially system-independent**: pp and the PbPb overlay agree to 0.13 % absolute
   (0.9133 vs 0.9121), and the unit-area pp/PbPb truth shape ratio is flat at 1 within the
   overlay's statistics across the whole axis. At truth level the two differ only through the
   beam mixture, and that difference is not resolvable here. (The overlay is a 108 k test
   sample and its panels are visibly noisy on the canonical fine axis; the binning was NOT
   coarsened for it.)

**Working point.** Deliberately absent — a truth muon has no reconstruction quality, so the
Tight/Medium config var every other plot set carries has no meaning. Recorded as such in
`docs/muon_wp_registry.md` §8 (new section: "WP-INDEPENDENT by construction").

### F13 — TLegend boxes were covering curve segments in the new-gap-cut figures (2026-08-13)

User-reported, on `muon_gap_cuts/new_gap_cuts/` (the legacy `old_gap_cuts/` set was explicitly
out of scope and was NOT regenerated, so it still carries the old layout). The legends are
opaque by necessity — they sit over the shaded gap bands — and are drawn last, so a legend at a
hard-coded corner erases whatever passes underneath it. The corner that is empty in pp is
crossed by the Pb+Pb acceptance edge, which is why the defect showed up on the Pb+Pb panels.

Three distinct faults, all fixed in `plot_muon_q_eta_spectrum.cxx`:
1. **Panel 1's legend carried the entire gap-window key** — five entries, the full width of the
   frame — so there was nowhere to put it that was not on top of the spectrum; it sat across the
   q·η ≈ 0 dip and the forward edge. The key now goes **once across the canvas** in a band freed
   above the pads (pads 0.90 → 0.875), and panel 1 has no in-frame legend at all: it has one
   curve, which the panel header already names.
2. **The pp-vs-Pb+Pb legend was sized for two entries but held five**, so ROOT printed the
   window lines on top of the sample lines. Same fix (key to the canvas band) plus a new
   `LegendHeight(n_entries, text_size)` so a box can no longer be sized independently of what
   it holds.
3. **The remaining 2–5 entry legends are now auto-placed** by `AutoLegendBox()`: it scores the
   four in-frame corners by how much of the drawn curve falls inside the candidate box and
   returns the emptiest. Scoring walks consecutive-bin **segments**, not bin contents, because a
   step histogram draws the vertical connector between bins — a steep plunge (exactly the
   forward acceptance edge) can cross a box without any single bin content landing inside it.

Both WPs regenerated (`new_gap_cuts/` + `new_gap_cuts/medium/`, 8 PNGs); every panel visually
re-checked. The same helpers were ported verbatim into
`plot_muon_truth_q_eta_spectrum.cxx` so the truth and data figures place their legends by the
same rule and can be read side by side.

### F14 — q·η shape vs single-muon p_T: the adopted windows are p_T-stable, but a p_T-dependent depletion at q·η ≈ +0.6 and +1.2 is NOT covered (2026-09-03)

**Step 8, user request.** New standalone macro
`plotting_codes/single_b_analysis/plot_muon_q_eta_pt_dependence.cxx` →
`muon_gap_cuts/new_gap_cuts/muon_q_eta_pt_dependence.png` (Medium WP would go to
`new_gap_cuts/medium/`; only Tight was produced). 6 panels, 3 rows × 2 cols: pp24, then
Pb+Pb 23+24+25 at 0-10 / 10-20 / 20-30 / 30-50 / 50-80 %. Each panel overlays 4 unit-area
q·η PDFs for p_T ∈ [4,4.5), [4.5,5), [5,6), [6,∞) GeV. Separate macro from
`plot_muon_q_eta_spectrum.cxx` so the four existing figures need no rerun.

**Binnings.** q·η = `ParamsSet::makeEtaTrigEffcyBinning(1)` (184 bins, [−2.40, 2.40]), the
same canonical axis as every other panel in this doc; centrality = `ParamsSet::ctrbins` with
0-5 + 5-10 merged, identical to the sibling's Plot C; windows read from
`ParamsSet::single_mu_fiducial_gap_cuts`. **The four p_T slices are the USER'S EXPLICIT
CHOICE for this figure** (2026-09-03) and are deliberately NOT
`ParamsSet::single_mu_pt_coarse_bins` ({4,8,14,25,100}), whose single 4-8 GeV bin cannot
resolve the behaviour just above threshold. Plot-local diagnostic axis; feeds no physics
result, no fit and no correction. Nothing retyped, nothing silently redefined.

**Result 1 — the p_T-INDEPENDENT fiducial cut is justified where its windows sit.** The
q·η ≈ −1.11 barrel/endcap dip, the η ≈ 0 crack and the forward collapse appear at the same
q·η and roughly the same depth in all four p_T slices, in pp and at every centrality. The
rejected fraction is correspondingly near-flat in p_T:

| panel | 4.0-4.5 | 4.5-5.0 | 5.0-6.0 | >6 GeV | all p_T |
|---|---|---|---|---|---|
| pp24 | 0.0368 | 0.0329 | 0.0301 | 0.0425 | **0.0379** |
| Pb+Pb 0-10% | 0.0578 | 0.0544 | 0.0524 | 0.0563 | 0.0555 |
| Pb+Pb 10-20% | 0.0586 | 0.0555 | 0.0534 | 0.0573 | 0.0565 |
| Pb+Pb 20-30% | 0.0574 | 0.0549 | 0.0542 | 0.0576 | 0.0563 |
| Pb+Pb 30-50% | 0.0570 | 0.0547 | 0.0549 | 0.0577 | 0.0564 |
| Pb+Pb 50-80% | 0.0575 | 0.0567 | 0.0555 | 0.0591 | 0.0576 |

Range 0.030-0.043 (pp) and 0.052-0.059 (Pb+Pb), no monotonic trend. The pp all-p_T value
**0.0379 reproduces the 3.79 % of F10/F11a exactly**, and there is no centrality dependence
worth speaking of (0.0555 → 0.0576 from 0-10 % to 50-80 %).

**Result 2 — OPEN CUT QUESTION. The p_T dependence lives on the POSITIVE q·η side, where the
adopted set has no window.** In every panel the 4.0-4.5 GeV curve carries deep, narrow
depletions at **q·η ≈ +0.6** and **q·η ≈ +1.1-1.3** that shrink monotonically through
4.5-5.0 and 5.0-6.0 GeV and are essentially absent above 6 GeV. Those are precisely the two
legacy `ParamsSet::charge_eta_gap_cuts` windows **{0.56, 0.67}** and **{1.064, 1.29}** — and
precisely why the legacy `PassSingleMuonGapCut` gated them at p_T < 6 GeV
(`RDFBasedHistFillingData.cxx:144`). So the legacy code's p_T gate was tracking something
real, it is still present in the data, and the adopted set (F10/F11a) does not cover it.
Physically this is the expected toroid bending-out of low-p_T muons of that sign at the
detector feet / support structures (KB `physics/detector/atlas_run2_muon_trigger.md` lists
"gaps at η≈0, feet, support structures"). **This is a cut decision, i.e. the user's**: adding
a fourth and fifth window at all p_T costs yield above 6 GeV where there is no dip, whereas
leaving them costs a p_T-dependent acceptance that ε_acc must then carry. NOT acted on.

**Result 3 — the pp forward region is trigger-shaped, not gap-shaped.** For q·η > +1.5 the pp
low-p_T slices fall by more than a decade toward the acceptance edge while p_T > 6 GeV stays
flat. No Pb+Pb counterpart. This is the **2mu4 both-leg** turn-on folded into the acceptance
edge (Pb+Pb is single `mu4`), not detector structure — consistent with F12, which finds the
truth spectrum featureless there.

**Plot-review outcome (partial).** `/review-plot` iteration 1 returned **FAIL** on one
CRITICAL: the opaque in-frame `TLegend` covered the η ≈ 0 crack minimum in **all six**
panels — **a recurrence of F13 in new code**, because this macro hard-coded the box instead
of porting `AutoLegendBox()`. Measured: legend top = 0.119 on the shared Pb+Pb axis / 0.0272
on the pp axis vs crack minima 0.058-0.084 (Pb+Pb, all four curves in every panel) and
0.0231 (pp, p_T > 6 GeV), so the crack read ≈ 0.44 of plateau where it is really ≈ 0.25 —
the figure understated the very feature the middle window is drawn around. **Fixed** by
deleting the in-frame legend and drawing both keys once across the top of the canvas
(`CanvasKey()`, `kPanelTop` 0.965 → 0.942); the defect is recorded in that helper's comment
so it is not reintroduced a third time. Also fixed: the trigger is now named on every panel
(`HLT_2mu4` / `HLT_mu4`, sub-line) — it had been on pp only; the sibling's |q·η| ∈ [1.4,2.2]
binning-step caveat added to the header (the apparent step at |q·η| ≈ 1.4 is the 0.02 → 0.10
step change, not structure); Medium WP routed to `new_gap_cuts/medium/`.

Iteration 1 independently re-derived every number through a **different read path**
(`TTree::Draw` expressions instead of `SetMakeClass` + branch addresses) and got MATCH on all
of them — windows, axis, centrality edges, the pp24 population 9 107 895, every per-slice and
all-p_T fraction above, `Integral("width") = 1.000000` for all 12 curves checked, and the
axis-cap bookkeeping (9 bins off scale in pp; the Pb+Pb cap verified inactive). That also
rules out the `SetMakeClass(1)` branch-type trap. C3 came back
`RUN2-CROSSCHECK UNVERIFIED` — no Run 2 analog of this figure exists; the KB corroborates the
gap structure qualitatively only. Non-blocking.

⚠ **CAVEAT: iteration 2 was stopped by the user before returning a verdict.** The amended
figure has therefore **not** been independently re-verified — the legend fix was checked by
the agent only (crack now visibly reaching ≈ 0.065 in Pb+Pb, ≈ 0.022 in pp, i.e. no longer
truncated). The numbers in Result 1 are unaffected: the amendment touched drawing code and
output paths only, not the fill loop, and those numbers were verified in iteration 1. If this
figure is ever promoted toward a note, re-run `/review-plot` on it first.

### F15 — q·η spectra (data + truth) switched to linear y; `old_gap_cuts/` deleted (2026-09-06)

**Step 9, user request.** Both q·η spectrum macros — `plot_muon_q_eta_spectrum.cxx` (data,
F6/F9/F11/F13) and `plot_muon_truth_q_eta_spectrum.cxx` (truth, F12) — switched from log-y to
LINEAR y (floor 0, ceiling 1.15× the overlaid peak via a new `SetCommonLinRange`, replacing the
log-floor `SetCommonLogRange`/`kMaxDecades`/`NoteOffScale` machinery), so the plateau and the
three gap-cut windows' depth/steepness read on an absolute scale instead of a compressed one.
`/review-plot` PASS at iteration 1 (log
`.claude/logs/review-plot-20260906-212055-qeta-linear-yscale.md`); every console number
(gap-cut fail fractions, truth ε_acc and its per-window breakdown) verified byte-identical to
F6/F9/F11a/F12 by the reviewer's own independent rerun. `muon_q_eta_pt_dependence.png` (F14) is
UNCHANGED — out of scope, still log-y, a separate diagnostic.

**Regression found and fixed during execution (data macro only).** Switching to linear-y made
ROOT invoke its "×10^n" axis-multiplier annotation (needed once raw counts render as full
numbers up to ~3.5M), which collided with the `Header()`/`CanvasGapKey()` text drawn just above
the frame — present on every panel with an empty-sub `Header()` call and on the sub-line of
two-line headers. Root-caused and fixed, not patched with fragile per-panel offsets:
`TGaxis::SetMaxDigits(7)` (full-digit labels, no multiplier box) + `StyleSpectrum`'s y-axis
`SetTitleOffset` 1.70→3.00 + `SetLeftMargin` 0.20→0.28 in `MakeSinglePad` and BOTH
`p_main`/`p_rat` of `MakeRatioPads` (kept equal so the stacked main+ratio x-axes still line up).
Verified clean on the highest-magnitude case (Medium WP PbPb combined, ~3.5M plateau) and the
tightest pad grid (Plot C, 2×3). The truth macro's values stay at O(10–50) nb (weighted
cross-section, not raw counts) and never trigger this, so it needed none of these changes —
log-y removal only.

**`old_gap_cuts/` deleted per user instruction** (superseded by F10's settled cut set; the
legacy `PassSingleMuonGapCut` comparison it existed to provide is no longer needed). The macro's
`use_new_gap_cuts=false` code path (which would regenerate it) is left in place, unused, in case
a future legacy comparison is wanted again.

All 19 regenerated PNGs (8 `old_gap_cuts/[medium/]` — since deleted, see above — + 8
`new_gap_cuts/[medium/]` + 3 `truth_q_eta/`) visually re-checked; no header/label/legend
collisions remain on any panel inspected.

### F16 — F14's p_T-dependence macro also switched to linear-y; F15's header fix redone with native scientific notation (2026-09-06)

**User follow-up to F15.** Two changes:

1. **`plot_muon_q_eta_pt_dependence.cxx` (F14) switched from log-y to linear-y**, matching the
   other two q·η macros — no longer out of scope. Its content is unit-area PDFs (max ~0.45), so
   it never triggers ROOT's axis-multiplier notation and needed no header/margin changes, only
   the standard log-y-removal treatment (`SetCommonLogRange`/`kMaxDecades`/`NoteOffScale` →
   `SetCommonLinRange`, floor 0, ceiling 1.15× peak). Fixed a latent bug found along the way:
   `DrawGapBands` unconditionally computed `std::pow(10., gPad->GetUymin()/GetUymax())` for the
   gap-band vertical extent — correct only for a log axis, silently wrong (though harmless while
   the axis WAS always log) for the now-linear one; replaced with the raw `GetUymin()`/
   `GetUymax()` values.
2. **F15's header/axis-multiplier collision fix redone**: the F15 approach
   (`TGaxis::SetMaxDigits(7)` forcing full-digit y-axis labels, e.g. "3000000", plus a
   `SetTitleOffset`/`SetLeftMargin` increase to fit them) is REPLACED per explicit user
   instruction — full-digit labels are not wanted; ROOT's native `×10^n` scientific notation
   (e.g. "3500" ticks with a "×10^3" box) is. Reverted `TGaxis::SetMaxDigits`/`TitleOffset`
   (1.70)/`LeftMargin` (0.20) to their pre-F15 values in `plot_muon_q_eta_spectrum.cxx`, and
   fixed the ACTUAL collision at its source instead: `Header()`'s near-frame text offsets raised
   (empty-sub 0.025→0.055; two-line sub-line 0.020→0.050, main 0.067→0.095) so the header text
   clears ROOT's native multiplier box rather than sitting on top of it.

`/review-plot` PASS at iteration 1
(`.claude/logs/review-plot-20260906-224247-qeta-ptdep-linear-and-header-fix.md`); reviewer
independently reran all three macros and matched every console number exactly (F11a/F14/F12).
Verified clean (native "×10^3" box visible, no collision) on the highest-magnitude case (Medium
WP PbPb combined, ~3.5M plateau) and the tightest pad grid (Plot C, 2×3), and on both
`plot_muon_q_eta_pt_dependence.cxx` PNGs (Tight + Medium).

### F17 — CUT SET CHANGED by the user; PAIR-LEVEL |eta^pair| < 2.2 added; full gap-cut inventory (2026-09-07)

**User instruction (three parts).**

1. `ParamsSet::single_mu_fiducial_gap_cuts` **{{-1.20,-1.05},{-0.10,+0.06},{2.30,2.40}} ->
   {{-1.30,-1.05},{-0.10,+0.06},{2.20,2.40}}**. The barrel/endcap window is widened on the
   negative-q·η side; the forward edge goes back to 2.20, where it again coincides with the
   historical standalone per-muon `q·η < 2.2` signal cut. The crack window is unchanged.
2. The q·η spectrum plot set is regenerated at the new windows, overwriting `new_gap_cuts/`
   and `truth_q_eta/`.
3. A **PAIR-LEVEL gap cut |eta^pair| < 2.2** (symmetric) is added, applied to the signal pairs
   and everywhere the single-muon windows are applied to BOTH legs. Physics reason, per the
   user: a pair can reach |eta^pair| > 2.2 with both muons passing the one-sided q·η windows,
   but only through a charge-dependent corner of phase space whose pair efficiency is shaped by
   the single-muon cut itself rather than by a smooth detector response. Removed, not modelled.

**FORCED consequence — the coarse turn-on binning moved with the forward edge.**
`RDFBasedHistFillingData::SetIOPaths` throws at startup unless
`CommonEffcyConfig::q_eta_proj_ranges_coarse_incl_gap`'s top edge EQUALS the forward window's
lower edge. Top bin **{2.0, 2.3} -> {2.0, 2.2}**. Not an independent binning choice — it is
slaved to the cut; without it every trigger-efficiency job aborts before reading an event.

**USER DECISION (asked, 2026-09-07):** the coarse **pair-η** binning
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap` outer bins **±(2.0, 2.4) -> ±(2.0,
2.2)**, so the 9 panels exactly tile the surviving region. Keeping ±2.4 would have left the two
outer panels populated over half their width — dσ/dη diluted ~2× against the interior panels,
with labels that still said (2.0, 2.4).

**USER DECISION (asked, 2026-09-07):** the pair-level cut is **NOT** applied in the data
tag-and-probe. ε^nc is a per-muon efficiency; the partner's position does not enter its
definition (which is also why the tag is left uncut). The pair cut lives only where a PAIR is
the object.

**New measured costs** (macros' own console output, this session):

| | pp24 | PbPb 23+24+25 |
|---|---|---|
| reco muons removed, Tight | **6.93 %** (was 3.79) | **9.03 %** (was 5.60) |
| reco muons removed, Medium | 7.03 % | 9.77 % |
| truth ε_acc | **0.8789** (was 0.9133) | **0.8765** (was 0.9121) |

Truth cost per window: `[-1.30,-1.05]` 5.52 %, `[-0.10,+0.06]` 4.16 %, `[+2.20,+2.40]` 2.44 %
(pp24 fullsim; PbPb overlay 5.58 / 4.15 / 2.62 %). PbPb is centrality-flat (9.00–9.15 % over
0–10…50–80 %). pp p_T dependence is mild and monotonic (5.30 % in 4.0–4.5 GeV -> 8.00 % above
6 GeV); PbPb flat. Widening the barrel/endcap window and pulling the forward edge back roughly
**doubles** the single-muon cost of the fiducial cut. **These are single-muon numbers only —
the pair-level |eta^pair| < 2.2 cut's own cost is NOT yet measured.**

**Regenerated plot sets** (overwritten, `plots/single_b_analysis/muon_gap_cuts/`):
`new_gap_cuts/muon_q_eta_spectrum_{pp24,pbpb_combined,pbpb_combined_ctr_binned,pp_vs_pbpb_shape}.png`
(+ `medium/`), `new_gap_cuts/muon_q_eta_pt_dependence.png` (Tight + Medium),
`truth_q_eta/` (3 PNGs).

**COMPLETE GAP-CUT INVENTORY** (two independent subagents, cross-checked; their scratch docs
`_sub_gapcut_inventory_{A,B}.md` merged here and deleted). The pre-existing summary in the
`ParamsSet.h` STATUS block was correct but **incomplete** — it missed 7 live sites.

*Sites reading `ParamsSet` (the fiducial cut), LIVE:*

| Site | Step | Legs | Pair cut added? |
|---|---|---|---|
| `PP.cxx` `FillHistogramsGeneric` | pp24 generic + MC-data compr | both | **yes** |
| `PP.cxx` `signal_cuts` | pp24 crossx signal region | both | **yes** |
| `PP.cxx` template-fit pass `signal_cuts_no_minv` | low-mass template | both | **yes** |
| `PP.cxx` template `_nosel` | low-mass template | both | **yes** |
| `PP.cxx` template inputs in the nominal crossx pass | low-mass template | both | **yes** |
| `PythiaFullsim.cxx` `gap_truth` / `gap_reco` | pp24 pair reco-eff | both, each leg | **yes** (truth_pair_eta / pair_eta) |
| `PowhegTruth.cxx` generic `_gapcut_truth` | POWHEG truth generic family | both (truth) | **yes** |
| `PowhegFullsim.cxx` `_single_b_pass_signal_truth_gapcut` | MC-vs-data signal family | both (truth) | **yes** |
| `FillMCTrigEffHists.cxx` `kGapLeg` | MC trig-eff Steps 2 & 4 | both legs | **yes** |
| `FillMCTrigEffHists.cxx` `kGapPair` | MC trig-eff Step 3 + kn-split stats | both | **yes** |
| `FillMCTrigEffHists.cxx` `kGapSingle` | MC trig-eff Step 1 + sanity | single muon | **no** — no pair exists |
| `Utilities/MCTrigEffPairSelection.h` `FiducialGapCut()` | closure (`FillMCTrigEffClosure.cxx`) + `SingleBSignalCutsReco()` | both | **yes** |
| `PP.cxx:172` / `PbPb.cxx:340` T&P probe | data mu4 ε^nc | probe only | **no** — user decision |

*Sites still on the RETIRED one-sided `q·η < 2.2`, numbers retyped — NO gap cut, so no pair cut
either (unchanged, and each is a pre-existing divergence, not one introduced here):*
`PbPb.cxx` crossx signal region + both template blocks; `PythiaFullsimOverlay.cxx`
per-centrality `pass_signal_{truth,reco}` (note the class is INTERNALLY SPLIT — its inclusive
nodes inherit the fiducial cut from the Pythia-fullsim base); `PythiaTruth.cxx` truth template
`_sigsel` and truth signal acceptance; `PowhegTruth.cxx` signal acceptance;
`PowhegFullsim.cxx` legacy `pass_signal_{truth,reco}`. Diagnostics likewise:
`plot_sig_accept_cutflow_above_60GeV.cxx`, `plot_dr_vs_pair_pt_diagnostic.cxx`,
`plot_reco_distr_singleb_vs_op_pp24.C`, `fill_pbpb_muon_pt45_diag_counts.cxx`.

*NTuple processing applies NO gap cut* (verified by grep across `NTupleProcessingCode/`); only
|η| < 2.4. Every consumer cuts for itself.

*Legacy `_wgapcut` family* (`PassSingleMuonGapCut`: |η| < `eta_gap_cut1` = 0.135 at all p_T +
`charge_eta_gap_cuts` below 6 GeV) is a DIFFERENT object, live only as a pp24 diagnostic
(`run_crossx_hist_filling_pp24.sh` sets `output_gapcut_hists=true`). Untouched.

**Open defects found by the inventory, NOT fixed here** (each is pre-existing):

1. `PbPb.cxx:707` loads the trigger-efficiency input measured WITH the probe fiducial cut
   (`..._qeta_fid...root`) while the PbPb signal selection is still `q·η < 2.2` — efficiency and
   selection on different fiducial regions.
2. R_AA therefore mixes fiducial regions (pp numerator fiducial, PbPb denominator not).
3. `ParamsSet::PassSingleMuFiducialGap` has **zero call sites**; the three q·η plotters
   re-implement the test with OPEN intervals in DOUBLE, whereas the helper and
   `FiducialGapCutExpr` use CLOSED intervals in FLOAT. Measure-zero, but it is exactly the bug
   class fixed at source on 2026-08-18.
4. Truth signal acceptance (Pythia + POWHEG) is on the retired cut, so it cannot yet supply
   ε_acc for the current fiducial region.
5. `PowhegFullSimExtras.c:35` `if ((muon.eta) > 2.4)` — missing `fabs` (obsolete path).

### F18 — the pair-η binning change silently truncated Pb+Pb; now it throws (2026-09-07)

Found by `/review-analysis-code` on the F17 change (CRITICAL, verified independently by bin
arithmetic). Shrinking the coarse pair-η panels to ±2.2 is safe on the **pp** axis
(`ParamsSet::N_PAIR_ETA_CROSSX_BINS = 48` over [−2.4, 2.4], width exactly 0.1 → 2.2 **is** a bin
edge) but **not** on Pb+Pb, which still books a retyped `44, −2.4, 2.4` axis (width 0.10909…)
and has adopted neither the fiducial nor the pair-level cut. On that axis
`FindBin(2.2 − 1e-6)` = bin 43, upper edge **2.29091**, where it used to be bin 44 (upper edge
2.4): every Pb+Pb pair with |η^pair| ∈ [2.29091, 2.4] would have disappeared from the panel set
and from the pair-η-integrated dσ/dη, with no warning, while the outer panel labelled (2.0, 2.2)
actually integrated 1.96364–2.29091.

**Fix applied — a guard, not a silent choice.** `Utilities/PairEtaPanelBins.h` (new):
`PairEtaPanels::Bins(axis, panel, context)` computes the bin range and THROWS unless both panel
edges are bin edges of *that* histogram's η axis, naming the range actually covered, the axis's
bin count/limits, and the fix. It also carries `CheckPanelsMatchFiducialCut()` — the pair-η
analogue of the q·η startup throw — which fires if the panel set ever stops tiling exactly
|η^pair| < `ParamsSet::pair_eta_fiducial_max` (that coupling was asserted in three comments and
enforced nowhere).

**The guard had to be HOISTED — this is the part the first fix got wrong.** Round 1 put it inside
`SingleBCrossxPlotterBase`, which covers five projection sites but not the *five more outside that
class* that project the same panels: `plot_crossx_trig_corr_sanity.C` (run by BOTH the pp and the
Pb+Pb crossx pipelines, reading the Pb+Pb 44-bin histograms), `plot_crossx_reco_eff_stages.C`
(same Pb+Pb family), `SignalAcceptancePlotter.cxx` (Pythia **and** POWHEG truth acceptance, both
booked `44, -2.4, 2.4`), `plot_pp_counts_pair_pt_in_eta.cxx` (a derived class that still called
`FindBin` directly), and `plot_mc_data_pair_pt_in_eta.cxx`. Found by the iteration-2 reviewer,
also CRITICAL. **A guard only some call sites use is not a guard** — all ten sites now go through
the header.

**A THIRD class of consumer, found in iteration 3 (also CRITICAL).** The two MC-trig-eff closure
plotters (`plot_mc_trig_eff_closure.cxx`, `plot_mc_trig_eff_closure_compare.cxx`) do not call
`FindBin` at all — they treat pair-η histogram **bin `iz`** as **panel `iz`** and label it
`EtaLabel(eta_ranges[iz-1])` read live from `CommonEffcyConfig`, guarded only by a **bin-COUNT**
check. A count check cannot see an edge move: all ten `mc_trig_eff_closure_*.root` files on disk
still carry 9 bins over ±2.4, so they would have passed and been drawn with panels labelled
(2.0, 2.2) that actually contain (2.0, 2.4). Both now call
`PairEtaPanels::CheckAxisAligned` per panel after the count check, so a stale file throws.
**Lesson: "project a bin range" was the wrong thing to grep for; "identify a bin with a panel" is
the actual hazard, and it can be spelled without `FindBin`.**

Also fixed in iteration 3: the MC-vs-data pair-η-**integrated** panel was labelled
`-2.4 < η^pair < 2.4`, read off the histogram axis, on a spectrum every input of which now carries
|η^pair| < 2.2 — a final-results figure overstating its own acceptance. It is now labelled from
`ParamsSet::pair_eta_fiducial_max` (the projection still runs over the full axis, which is
correct: nothing is left out, the region is simply empty).

`plot_pp_counts_pair_pt_in_eta.cxx`'s own partition assertion was also repaired: it asserted that
the nine panels tile the whole ±2.4 fine axis, which they no longer do, and its throw message
misdiagnosed the failure as `N_PAIR_ETA_CROSSX_BINS` being wrong. It now says what a residual
actually means — an input file filled before the pair-level cut existed.

**Unblocking Pb+Pb is a user decision with a rerun:** move `RDFBasedHistFillingPbPb.cxx` onto
`ParamsSet::N_PAIR_ETA_CROSSX_BINS` **and** onto the fiducial + pair-level cuts, then rerun the
Pb+Pb crossx hist filling — in one step, or its code and its on-disk histograms disagree. The same
applies to the Pythia/POWHEG **truth signal-acceptance** producers, which also book the 44-bin
axis; until then the signal-acceptance panel plots throw too.

**Iteration 4 — the hazard class is CLOSED, and two guards were hardened.** The reviewer's own
three-pass sweep (greps for `FindBin`, `Projection[XYZ]`, `GetBinContent` loops, suffix read-back,
`RangesToEdges`, and bin-index-as-panel-index, across all of `plotting_codes/`, `Utilities/`,
`RDFBasedHistFilling/`, `pipelines/` and the loose top-level macros) found **no remaining
unguarded panel→bin association**. Every site is either guarded by `PairEtaPanelBins.h`, carries
its own throwing whole-axis edge comparison (`PairRecoEffEvaluator`, both
`DrCorrection*Evaluator`s, `plot_pp24_fullsim_pair_reco_eff`, `FillMCTrigEffClosure`,
`fit_dr_corrections`), merely BOOKS an axis from the ranges, or is structurally immune (bins and
labels both read from the file's own axis). Two guards were still too weak and are now fixed:

* **The MC-vs-data integrated panel asserted the fiducial region without checking it.** Its label
  now says |η^pair| < 2.2 while the projection deliberately runs over the whole ±2.4 axis, so the
  claim is only true if the inputs carry the cut — and the reviewer measured **1.66 %** of the
  current pp24 input outside it. `assert_inside_fiducial` now throws on each of data / Pythia /
  POWHEG. `PairEtaPanels::Bins` guards ALIGNMENT; this guards EMPTINESS — they are different
  properties and the figure needs both.
* **`build_pp24_fullsim_pair_reco_eff.C` was the one WRITER in the chain with no guard.** It
  stamps a provenance string claiming the panel axis and the pair-level window; the current
  `pair_reco_eff_pp24_full.root` on disk has the OLD ±2.4 outer bins while already carrying that
  claim. Both consumers throw, but only after a mislabelled file exists. The writer now checks
  its own axis before writing.

**Also found in iteration 4, PRE-EXISTING and now made loud, not fixed:** `RAA_plotting.cxx`
mode 2 (R_AA vs pair η) forms an **index-wise** PbPb/pp ratio while Pb+Pb has 44 pair-η bins and
pp has 48 — at bin 44 that divides Pb+Pb [2.291, 2.4] by pp [1.9, 2.0], and every bin above the
first is offset. The in-code comment claimed the axes "match by construction", which is true for
mode 1 (both `pT_bins_120`) and false for mode 2. This predates the 2026-09-07 change; its real
fix is the same Pb+Pb migration to `N_PAIR_ETA_CROSSX_BINS` + rerun. A bin-count and per-bin
edge check now **throws** instead, and the display range is taken from
`ParamsSet::pair_eta_fiducial_max` rather than a retyped 2.3.

**One more doc defect the same review caught:** `signal_selection_change_impact.md` §1 listed the
trigger-efficiency derivation as UNCHANGED — true for a *pair*-cut change, but this change moved
the **probe** gap windows and the coarse q·η binning, so P2/P3 MUST be refit. §1 now carries an
explicit carve-out.

Also fixed in the same amendment round (all WARNING-level, all comment/provenance truthfulness):
the `ParamsSet.h` header still called the cut "PROVISIONAL / not yet implemented / LIVE in the
trigger-efficiency chain only" and still quoted ε_acc = 0.9133 twenty lines below the new table;
`plot_muon_q_eta_spectrum.cxx` carried a literal window list stale by two revisions (now removed
in favour of the run-time printout); `plot_forward_qeta_edge_scan.cxx` still said the edge "was
loosened" to (2.30, 2.40); `FillMCTrigEffHists.cxx` labelled the pair-η axis "[-2.4,2.4]";
`dr_correction_cell_groups.h` AND `pipelines/run_dr_correction_fits.sh` documented the etamerge
top group as "2.0 ≤ |η| < 2.4"; and three
provenance strings **written into produced artefacts** (the Step-3 pair-statistics CSV, the
ε_reco histogram title, the raw-counts CSV) described only the single-muon windows — all three
now also emit `ParamsSet::PairFiducialEtaCutExpr`. `FillMCTrigEffHists.cxx`'s `kGapPair` was
migrated onto `Utilities/MCTrigEffPairSelection.h::FiducialGapCut()`, honouring that header's
MIRROR NOTICE, so the gap construction can no longer drift between the two.

### F19 — SUPERSEDED: ε_acc is no longer a separate factor (2026-09-17)

User decision, implemented in `pair_reco_eff_gap_acceptance.md`: the gap cuts (single-muon q·η
windows AND the pair-level |η^pair| < 2.2) are DETECTOR-ACCEPTANCE cuts, not signal cuts, and
their loss is carried by the pair reconstruction efficiency ε_reco — gap cuts on the RECO leg
(numerator) only, truth denominator = signal cuts alone. The "ε_acc itself — not built" item of
this doc's Latest Stage is therefore closed as superseded, not built: the F12/F17 truth
ε_acc numbers (0.9133 / 0.8789) are historical and must NOT be applied — they are already inside
ε_reco, and applying them would double count. The site inventory (F17) is unchanged; only the
pp24 fullsim `pass_signal_truth` row is now "reco leg only".

## Ruled Out (append-only)

- *Plotting from the existing 2026-06-23 trees* — stale on the one-sided Δp/p fix (F3), which
  is η-dependent and therefore not safely ignorable for a gap-region study.
- *Re-deriving the spectrum standalone from the raw NTUPs* — forbidden by the NTuple-Processing
  Provenance rule and unnecessary: the single-muon trees are exactly the intended output.

## Latest Stage

**STEP 10 IN PROGRESS (2026-09-07, F17) — CUT SET CHANGED + PAIR-LEVEL CUT WIRED IN.**
`single_mu_fiducial_gap_cuts` = **{{-1.30,-1.05}, {-0.10,+0.06}, {2.20,2.40}}**, plus a NEW
**pair-level |eta^pair| < 2.2** (`ParamsSet::pair_eta_fiducial_max`), wired into every site
where the single-muon windows are applied to BOTH legs (12 sites; full inventory in F17).
Two coupled binning edges moved with it: the coarse q·η turn-on top edge 2.3 -> 2.2 (FORCED by
the startup throw) and the coarse pair-η outer bins ±(2.0,2.4) -> ±(2.0,2.2) (user decision).
q·η plot set regenerated; new costs 6.93 % (pp) / 9.03 % (PbPb) of reco muons, ε_acc = 0.8789.

**NOT YET DONE — the whole downstream chain is STALE.** Nothing has been rerun. Per
`docs/signal_selection_change_impact.md` this needs, in order: (1) the single-muon mu4 turn-on
refits (the coarse q·η binning changed AND the probe gap cut changed), (2) MC trig-eff Steps 1-4
+ the ΔR correction fits + closure, (3) pp24 fullsim reco-eff + rebuild of
`pair_reco_eff_pp24_full.root`, (4) pp24 data crossx + template inputs, (5) POWHEG truth/fullsim
generic+signal families, (6) all crossx / R_AA / MC-data-comparison / acceptance plots, (7) note
figure re-sync. **Awaiting the user's go-ahead on how much of that to run and in what order.**

**STEP 9 DONE (2026-09-06, F15+F16).** q·η spectra (data + truth + F14's p_T-dependence
diagnostic) all switched to linear y-scale, `/review-plot` PASS on both rounds, all numbers
verified byte-identical to F6/F9/F11a/F12/F14 by the reviewer's own independent reruns.
`muon_gap_cuts/old_gap_cuts/` deleted per user instruction (superseded by F10's settled cut
set). The header/axis-multiplier-box collision the linear switch exposed is fixed by
repositioning `Header()`'s near-frame text (F16) so ROOT's native `×10^n` scientific-notation
axis labels display normally — NOT by forcing full-digit labels (F15's approach, reverted).

**STEP 8 DONE (2026-09-03).** The p_T-dependence figure exists at
`muon_gap_cuts/new_gap_cuts/muon_q_eta_pt_dependence.png` (F14). Headline: the three adopted
windows are **p_T-stable** — the rejected fraction varies only over 0.030-0.043 (pp) /
0.052-0.059 (Pb+Pb) across the four p_T slices, with no trend — so the p_T-INDEPENDENT
fiducial cut is justified *where its windows sit*. Its plot review is only **partially**
complete: iteration 1 caught and fixed a CRITICAL recurrence of F13 (legend covering the
crack in all six panels) and verified every number through an independent read path;
iteration 2 was stopped before a verdict, so the amended figure is not independently
re-verified.

**NEW OPEN CUT QUESTION (F14 Result 2) — for the user.** The p_T dependence is on the
POSITIVE q·η side, where the adopted set has no window: the 4.0-4.5 GeV slice is deeply
depleted at **q·η ≈ +0.6** and **≈ +1.1-1.3**, shrinking monotonically with p_T and gone
above 6 GeV. Those are exactly the legacy `charge_eta_gap_cuts` windows {0.56, 0.67} and
{1.064, 1.29}, which the legacy code applied only below 6 GeV — so that p_T gate was
tracking something real and the adopted set does not cover it. Adding windows there costs
yield above 6 GeV where there is no dip; leaving them makes the acceptance p_T-dependent and
pushes the cost into ε_acc. **Not acted on — a cut decision is the user's.**

**STEP 6 + 7 DONE (2026-08-13).** The TRUTH q·η reference spectrum exists at
`muon_gap_cuts/truth_q_eta/` (F12) and the legend-over-curve defect in the new-gap-cut data
figures is fixed (F13). Headline: the −1.11 dip, the η ≈ 0 crack and the one-sided forward
collapse are **100 % detector** — the truth spectrum is featureless across all three
(1.036 / 1.007 / 1.005 where the data sit at 0.32 / deep / 0.094) — and the truth-level
acceptance is **ε_acc = 0.9133 (pp) / 0.9121 (PbPb overlay)**, i.e. the cut costs **8.67 %** of
TRUE muons against only 3.79 % of reconstructed ones. Building ε_acc itself is still the
deferred item below; F12 supplies its numerator/denominator.

**Cut set FIXED (F10, F11a):** `{{-1.20,-1.05}, {-0.10,+0.06}, {2.30,2.40}}` in q·η, at a cost
of 3.79 % (pp) / 5.60 % (PbPb) of single muons. Plots regenerated at the new window with the
obsolete no-fit dashed markers removed (F11c).

**IMMEDIATE, OPEN — the trigger efficiencies are STALE (F11b).** Unlike every earlier edit to
this vector, the 2026-08-12 crack widening DOES move results: `single_mu_fiducial_gap_cuts` has
been live in the trigger-efficiency chain since 2026-08-04 (MC trig-eff Steps 1–4; data
tag-and-probe probe leg, pp + PbPb). The MC and data trigger efficiencies must be refit.
**Not rerun here** — that pipeline belongs to the sibling session's `mc_trigger_efficiency.md`
and a unilateral rerun would collide with its in-flight work. Awaiting the user's call on who
runs it.

**Still deferred, awaiting an explicit go-ahead:**

1. **Wire the vector into the SIGNAL selection** — replaces the standalone per-muon
   `q·η < 2.2` at the 13 call sites (F8); full `signal_selection_change_impact.md` blast
   radius. Prerequisite from the 2.30 edge (F10.1) is already satisfied: the coarse binning
   now runs to 2.3.
2. **Reconstruction efficiency, data ntuple processing and the template fits** — the other
   three consumers named in the `ParamsSet.h` STATUS block, none of which apply the cut yet.
3. **ε_acc itself** — SUPERSEDED 2026-09-17 (F19): folded into ε_reco, never a separate factor.

Also open from F8: `q_eta_proj_ranges_fine_excl_gap` holes still do not align with the chosen
windows — but that binning is no longer what the trigger-efficiency fits use, so it is now a
tidy-up, not a correctness issue.

Side flags, never part of this task (F4): cross-year FCal scaling is documented but never
implemented in code; muon-level `ev_centrality` bypasses the PbPb25 pair-level centrality
recalculation.
