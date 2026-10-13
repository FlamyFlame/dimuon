# Analysis Overview — Run 3 Single-b Dimuon

**Purpose & scope of this document.** This is the high-level, *stable*
conceptual reference for the analysis: what we measure, why, and the
**scientific methodology**. It is meant to guide both implementation/
investigation sessions and academic writing (IntNote / paper). It deliberately
carries **no status, dates, or task tracking** — those live in the roadmap. If a
physics statement here ever conflicts with code or another doc, this document
(plus the relevant Physics Procedure in a tracking doc) is the intended
ground truth; flag the conflict rather than silently following the code.

**See also:** this document owns the *what/why*; the roadmap
(`docs/tracking/analysis_roadmap_2026_06.md`) owns *status / next steps*; code
architecture is in `README.md`. Full doc catalog: project `CLAUDE.md`
(Documentation References).

> **Terminology.** "Methodology" in this document means **scientific
> methodology** — the physics-analysis method. The three-stage code pattern
> (NTuple processing → RDataFrame histogram filling → plotting) is an
> *implementation design*, not the methodology; see `README.md`.

---

## 1. Objective & physics motivation

This is an ATLAS **Run 3 single-b → dimuon** measurement. We select muon pairs
produced in the semileptonic decay of a **single b-hadron** and use them as a
**proxy for the parent b-hadron** kinematics. By comparing Pb+Pb to pp we probe
**b-quark in-medium energy loss** (and, potentially, b-quark **azimuthal
anisotropy / flow**) in the quark–gluon plasma — i.e. how the heavy b quark is
modified by the medium relative to a pp baseline.

Datasets (Run 3): **pp 2024** (√s = 5.36 TeV) as the reference, and **Pb+Pb
2023 / 2024 / 2025 / 2026** (√s_NN = 5.36 TeV), combined for all Pb+Pb results
(Σ L_int = 6.98170 nb⁻¹; 2026 = 2.35549 nb⁻¹ = updated-GRL lumitable minus the D8 runs 522200/522041, 2026-10-13). Run 2 code paths (pp17, PbPb15/18)
exist only for cross-checks and are not maintained for new analysis decisions.

## 2. The signal: a single-b dimuon

The signal muon pair comes from **one b-hadron's decay chain**: one muon from the
direct semileptonic decay (b → μ X) and the other from the sequential charm decay
in the same chain (b → c → μ X). Because both muons share a single b parent, the
pair has characteristically **low invariant mass** and **small opening angle**.

**Truth-level signal region** (as implemented):
- opposite-sign (OS) muon pair,
- **both muons p_T > 4.5 GeV** (4 -> 4.5 on 2026-09-08). Applied UPSTREAM, at the NTuple stage,
  so every consumer inherits it; it is mode-dependent in data, the trigger-efficiency mode keeping
  4.0 GeV because ε^nc is a per-muon efficiency evaluated only above 4.5 (decision D8).
- m_μμ ∈ **[1.08, 2.9] GeV** (`ParamsSet::signal_minv_min/max` + `SignalMinvCutExpr()`, the
  SINGLE SOURCE since 2026-09-17 — every live data / MC-analog / trigger-efficiency site reads it).
  **Coupled consumers:** the pp24 single-value pair trigger efficiency applied above 74 GeV is
  MEASURED inside this window (`docs/tracking/pp24_trig_eff_hybrid_application.md` D3), and the
  low-mass template fit is to be performed inside this SAME window (no φ / J/ψ peaks; user
  decision 2026-09-17, not yet implemented). Changing the window ⇒ re-measure the pair
  efficiency (`run_mc_trigeff_pair_eff.sh`) and refill every consumer,
  `docs/signal_selection_change_impact.md`.
- pair p_T > **9 GeV** (`ParamsSet::signal_pair_pt_min`; 8 -> 9 on 2026-09-08),
- **Detector-gap fiducial cut**, `ParamsSet::single_mu_fiducial_gap_cuts` (2026-09-07):
  BOTH muons must have q·η outside {(-1.25,-1.05), (-0.10,+0.06), (2.20,2.40)} (closed
  intervals, one-sided because the windows track the toroid bending direction), **and**
  the pair must satisfy **|η^pair| < 2.2** (symmetric, pair-level). This replaced the
  former one-sided per-muon `q·η < 2.2` on 2026-08-17 in pp24; **Pb+Pb was brought over on
  2026-09-08**, so pp and Pb+Pb share the signal region again (the Pb+Pb histograms must be
  refilled before R_AA is quoted). **The gap cuts are DETECTOR-ACCEPTANCE cuts, not signal
  cuts** (user decision 2026-09-17, `docs/tracking/pair_reco_eff_gap_acceptance.md`): they are
  applied to the RECONSTRUCTED pair only, and their loss is carried by the pair reconstruction
  efficiency ε_reco^pair (§4b), whose truth denominator carries the signal cuts alone. The
  corrected cross-section therefore refers to the truth region WITHOUT single-muon gap
  windows: muon p_T > 4.5 GeV, |η| < 2.4; mass window; pair p_T > 9 GeV; and
  **|η^pair| < 2.2 as the measured range** — the coarse η^pair axis on which ε_reco is measured
  tiles exactly [−2.2, 2.2], so truth pairs outside it enter no efficiency cell (3.35 % relative
  to the in-range truth signal-region weight, 3.2 % of the whole region; their reconstruction
  inside the range is 0.02 % migration). The
  pair-level 2.2 boundary is thus a truth fiducial edge of the measurement, while the single-muon
  q·η windows are acceptance-corrected. There is **no separate acceptance factor ε_acc** (the historical ε_acc numbers in `muon_gap_cuts_acceptance.md` /
  `ParamsSet.h` are superseded and must not be applied).
- **No ΔR cut.** *(A former ΔR(μ,μ) > 0.05 cut was removed 2026-06-22 — interim
  nominal. It was a stats workaround for the data-based dR trigger-efficiency
  inverse-weighting, not a physics requirement, and it biased the high-pair-pT
  signal acceptance, α dropping ~0.72→0.38 by p_T≈100 GeV as the muons collimate
  below ΔR=0.05; removing it flattens α to ~0.72. The dR trigger correction is
  moving to MC; whether any ΔR cut is needed is decided once that MC exists. See
  `docs/tracking/remove_dr_cut_signal_selection.md`.)*

This single-b topology is what distinguishes the signal from the dominant
heavy-flavour background where the two muons come from **two different** heavy
quarks (see §6).

> **Changing this selection?** Any add/remove/value change to a cut here (incl.
> selection systematics) has a wide rerun blast radius. **Before editing, read
> `docs/signal_selection_change_impact.md`** (the full recompile→rerun→replot map).

## 3. Final observables & physics questions

1. **Differential production cross-sections** of single-b dimuons, in pp24 and in
   Pb+Pb (years combined), differential in **(pair p_T, pair η, centrality)**.
2. **Nuclear modification factor R_AA** — the primary result — quantifying b-quark
   suppression in the medium:

   R_AA = (1 / ⟨T_AA⟩) · (dN_PbPb / dX) / (dσ_pp / dX),  X ∈ {pair p_T, pair η}

   per centrality class, where ⟨T_AA⟩ is the Glauber nuclear overlap.
3. **Azimuthal anisotropy (e.g. v₂)** of the dimuons — *potential* additional
   observable, as a proxy for b-hadron flow. **In scope is an open question**
   (also whether a Run-2-style Δφ away-side correlation-width measurement is
   included); see roadmap §Q2.9. The flow-fit machinery is not assumed here.

**Physics questions:** Does the b quark lose energy in the QGP, and how does the
loss depend on p_T, η, and collision centrality? Does the b quark flow with the
medium?

## 4. Scientific methodology (physics method)

The measurement proceeds as a corrected, normalized yield:

```
dσ/dX (or dN/dX) = [ raw OS yield in signal region, background-subtracted ]
                   × (per-pair efficiency correction)
                   / (normalization)
```

**(a) Yield & combinatorial background.** Count OS pairs in the signal region.
The uncorrelated combinatorial component is estimated from the **same-sign (SS)**
pairs of the same sample (sign1 = SS, sign2 = OS in the trees). Residual
correlated backgrounds are handled by the purity/template procedure in §6.

**(b) Efficiency correction (per pair).** Each surviving pair is weighted by the
inverse of the product of two efficiencies:

```
w⁻¹ = ε_trig^pair · ε_reco^pair(pair p_T, pair η, ΔR)
```

- **ε_trig^pair** — *trigger* efficiency, built from the single-muon
  no-correlation efficiency ε^nc(p_T, q·η) measured in data:
  - PbPb single-mu4: inclusion–exclusion, ε_trig^pair = ε₁^nc + ε₂^nc − ε₁^nc·ε₂^nc;
  - pp24 2mu4 (HYBRID since 2026-09-17, `docs/tracking/pp24_trig_eff_hybrid_application.md`
    §2; `Utilities/PairTrigEffCrossxEvaluator.h`), over the canonical coarse pair-pT axis
    `ParamsSet::pair_pt_coarse_bins`:
    - bins 1–6 ([9, 74.24) GeV): ε_trig^pair = ε₁^nc · ε₂^nc · ε_ΔR(ΔR; cell), ε_ΔR the MC
      2mu4 close-by correlation correction from the pp24 fullsim (opposite sign, 3-group
      |η^pair| fold, free-baseline fit; expo, polynomial in three named forward cells,
      interpolation fallback; = 1 for ΔR ≥ 1);
    - bins 7–8 ([74.24, 150) GeV): ε_trig^pair = ε^pair_MC(cell) · SF₁ · SF₂ — the SINGLE-VALUE
      MC pair 2mu4 efficiency per (pair-pT, |η^pair|) cell in the signal mass window, a
      PAIR-level efficiency (NOT multiplied by ε^nc), corrected to data by the product of the
      two single-muon data/MC scale factors SF = ε^nc_data/ε_MC evaluated per leg.
    TEMPORARY (awaiting requested high-pT / same-sign-filtered MC): the refused 3-pair cell
    [105.5,150)×|η| [2.0,2.2) uses the pT-merged cell of its |η| group; same-sign pairs are
    weighted with the opposite-sign numbers.
- **ε_reco^pair** — *reconstruction* efficiency, a single **pair** efficiency
  binned in **(pair p_T, pair η, ΔR)**, measured in **Pythia fullsim** (pp) and
  **Pythia fullsim HIJING overlay** (PbPb, per centrality). It is **not** the
  product of two single-muon efficiencies: the two close-by muons are correlated
  at reconstruction (shared ID hits, ambiguity resolution), and the ΔR axis
  captures that correlation. In central Pb+Pb the efficiency is genuinely lower
  and centrality-dependent (real occupancy effect from the underlying event), not
  a bookkeeping artifact.
  **Superseded by the 3D bin-by-bin correction of §4c (2026-10-09; not yet in the crossx):** the
  applied cells become (pair p_T, η^pair, m_μμ) in RECO kinematics; ΔR is dropped as an axis.
  **Definition (2026-09-17; pair category 2026-10-06):** ε_reco^pair(cell) = N[truth
  opposite-sign pair (any origin) in the truth signal region, both muons reconstructed, WP, and
  the reconstructed pair passing the full data selection *including the detector-gap cuts*] /
  N[truth opposite-sign pair in the truth signal region — signal cuts only, **no gap cut**].
  All truth OS pairs, not only `from_same_b`: ε_reco is a detector property of the pair
  kinematics, and the origin test found no cell-level difference between single-b and other
  OS pairs in populated cells (+18 % MC statistics;
  `docs/tracking/pair_reco_eff_origin_test_and_views.md` R1). It thus **includes the gap-cut acceptance**: a
  truth pair reconstructed into a gap window or across |η^pair| = 2.2 is an inefficiency. A
  truth pair outside the measured η^pair range (|η^pair| ≥ 2.2 — the truth fiducial edge of
  the measurement, see §2) reconstructed inside it is bin migration for the detector-response
  step, not part of ε_reco
  (`docs/tracking/pair_reco_eff_gap_acceptance.md`).

**Signal acceptance A_sig (2026-10-05/06, `docs/tracking/signal_acceptance_mass_window.md`).**
The result is a detector-corrected (3D bin-by-bin, §4c) cross-section in truth pair p_T > 9 GeV, |truth η^pair| < 2.2 — these
two cuts DEFINE the measured range and are never acceptance-corrected. TARGET (decided 2026-10-06):
A_sig(pair p_T, η^pair) = N[truth single-b OS in the measured range, both muons p_T > 4.5 GeV &
|η| < 2.4, **truth mass window**] / N[truth single-b OS in the measured range, both muons |η| < 2.5
(ID coverage), **no muon-p_T cut, all masses**], applied LAST to the detector-corrected yield. It contains NO
gap cut (that acceptance is in ε_reco^pair above). The single-muon part needs a NEW private Pythia
truth sample without a muon-p_T generator filter (central samples filter at 3.7 GeV, the old
private one at 3.5 GeV) — NOT yet produced; currently only the mass-window factor exists (interim):
Pythia integrated A_sig^mass = 0.749 (≈0.76 → 0.71 with pair p_T, flat in η^pair); POWHEG (bb, no-resonance-cut NTuple) 0.513 — the difference is the B → J/ψ X → μμ fraction of single-b pairs (31 % POWHEG vs 3 % Pythia), an OPEN generator/signal-definition question (tracking doc R4).

**(c) Detector response — 3D BIN-BY-BIN correction (decided 2026-10-09,
`docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md`).** Bin migration is corrected
**bin-by-bin, not by a full (Bayesian) unfolding**, because the migrations are small: in the pp24
Pythia fullsim every 3D cell with N_eff ≥ 100 has bin-by-bin efficiency ≥ 0.81 and purity ≥ 0.82
(N[T]-weighted means 0.92), and the response is near-diagonal
(`docs/tracking/pair_bin_migration_bin_by_bin_feasibility.md` R1–R3). The 3D variables are
**(pair p_T, η^pair, m_μμ)** — **m_μμ replaces ΔR** as the third axis (285/288 cells populated vs
144/288 for ΔR, which is kinematically tied to p_T via ΔR ≳ 2m/p_T; ε varies smoothly, ≈4–5 %,
across the 4 uniform mass bins; out-of-range feed-in ≈2× smaller). The reconstruction-efficiency and
migration corrections are then ONE factor per RECO cell j:

```
N^T_j = Σ_{data pairs in reco cell j} w_trig / ε^reco_j ,   ε^reco_j = N[R_j] / N[T_j] = 1/C_j
```

with T_j = truth OS pairs in truth cell j passing the truth signal selections (no gap cut) and R_j =
truth-matched pairs passing WP + the full data reco selection (incl. gap cuts) in RECO cell j, truth
ANYWHERE — the feed-in from outside T is present in data, so it must be in the numerator for the
output to be N[T_j] (user decision 2026-10-11). The output (summed over m) is the
efficiency-corrected truth (p_T, η^pair) distribution | truth signal selections — exactly the
numerator population of A_sig, which is applied next. **Propagated to the pp24 cross-section, the MC–data
comparison and R_AA on 2026-10-13** (`docs/tracking/crossx_raa_reco_binned_propagation_20261013.md`; it replaces
the truth-binned ε_reco(p_T, η^pair, ΔR) of (b) there); remaining open: prior (truth-shape) dependence of C_j; feed-in from truth muons below the NTuple truth
fiducial (pT < 4.5 / |η| > 2.4) is not in the truth-seeded MC tree. The ε^reco η^pair axis stays the canonical
coarse one even though the single-muon gap windows carve narrow holes into η^pair inside some coarse
cells (user decision 2026-10-13, option A): the coarse-cell yield is correct only as far as MC models
the within-cell reco η^pair shape, so a **within-cell η^pair shape systematic** is assigned — reweight
the MC truth η^pair shape inside each (p_T, η^pair) cell so MC reco matches the trigger-corrected data
in 0.1 bins and recompute the cell ε (|δ| up to 5 % in the central cell below 26 GeV, ≤ 0.5 % in cells
without a gap window; tracking doc §3d/R4).

**Data/MC scale factors (2026-10-13).** ε^reco_j is an MC efficiency; it is brought to data with the
single-muon data/MC efficiency SFs of the skim's MuonEfficiencyScaleFactors tool (WP-matched), as a
pair = product of the two legs, applied on the MC reco numerator:
ε^{reco,SF}_j = Σ_{R_j} w·SF₁·SF₂ / N[T_j] (Tight: ⟨SF₁SF₂⟩ = 0.931, inclusive ε 0.592 → 0.551).
The tool gives no SF below p_T = 5 GeV (all muons in [4.5, 5) GeV); **user decision 2026-10-13: such a muon
takes the SF of the nearest valid slice, [5, 6) GeV, at its η** (MC-weighted mean); SF = 1 there is the
systematic variation (−3.1 % on ε at pair p_T 9–12.8 GeV, −1.0 % at 12.8–18.2, ≲0.5 % above; cancels in R_AA).

**Pb+Pb reco efficiency (placeholder until the HIJING-overlay ε exists; 2026-10-13).** The Run-2 Pb+Pb
single-muon MC fits (ε₁ε₂ per centrality and q·η) know neither the Run-3 gap-cut acceptance (folded into
ε_reco) nor the SFs, so their η^pair dependence is wrong for our signal region. **User decision 2026-10-13
("anchor to pp24"):** ε_PbPb(j, c) = C_Run2(p, c) · ε^{reco,SF}_j with C_Run2(p, c) = E_Run2(p, c) /
E_Run2(p, 60–80 %), E_Run2 = the Run-2 ε₁ε₂ averaged over the pp24 fullsim truth signal pairs of pair-p_T
bin p. Run 2 supplies ONLY the centrality dependence relative to 60–80 %; level, p_T, η^pair (gap
acceptance), m dependence and SFs are the pp24 correction's, so no Run-2-vs-Run-3 efficiency-level factor
enters R_AA (the per-centrality R_AA efficiency ratio is exactly 1/C_Run2). The Pb+Pb trigger correction
stays the placeholder single-muon union without a correlation term.

**(d) Normalization.** pp: integrated luminosity. Pb+Pb: number of events and the
Glauber ⟨T_AA⟩ per centrality (and the total hadronic cross-section for centrality
calibration). These external numbers must be confirmed before entering results
(roadmap §Q2).

**(e) Centrality (Pb+Pb).** From forward calorimeter (FCal) ΣE_T with 2023 Glauber
thresholds for every Run-3 year (2024 carries its own provisional threshold vector,
which is in fact identical to the 2023 one). For 2025 and 2026 the centrality percentile is recomputed
from FCal E_T at pair level, which **enforces** the 2023 calibration regardless of what the
skim wrote (registered user decision, 2026-09-10). This is an override, not a repair: the
2025 `centrality` branch was measured on 2026-09-11 to be populated and already on the 2023
calibration (mean 20.454, max 84, 4.08 % zeros; agreeing with the 2023-table recompute in
299 992 of 300 000 events), correcting an earlier claim that it was unfilled. **No
cross-year FCal scaling is implemented** — `fcal_scale_pbpb_20YY.root` appears in
documentation only, never in a source file (see `placeholder.md`).

**(f) R_AA & systematics.** Form R_AA from the corrected, normalized pp and PbPb
spectra; propagate efficiency, bin-by-bin detector-correction, normalization, and selection
systematics.

> The code realizes (a)–(f) through the three-stage pipeline (NTuple processing →
> RDF histogram filling → plotting). That staging is an implementation choice and
> is documented in `README.md`; it is not part of the methodology above.

## 5. MC samples and their roles

| Sample | Role |
|--------|------|
| **Pythia truth** (5.36 / 5.02 TeV) | Generator-level signal/flavour & origin studies |
| **Pythia fullsim pp24** | **Reconstruction efficiency** and **detector response** (3D bin-by-bin correction, §4c) (pp) |
| **Pythia fullsim HIJING overlay** | **Reconstruction efficiency** (PbPb), per centrality; overlay puts HIJING into the detector at hit level |
| **Powheg truth (bb / cc)** | **NLO signal template** for the template fit (see §6). **bb and cc are two generator modes and their cross sections are ALWAYS ADDED** — see below |

**POWHEG is always bb + cc.** The two are POWHEG *generator modes*, each requiring a b-b̄ (resp.
c-c̄) pair in the NLO hard scattering, so they are distinct non-overlapping contributions to one
process and σ = σ_bb + σ_cc everywhere POWHEG appears (flavour/origin plots, MC-vs-data, the NLO
template). It holds even for observables that select only b-decay muons: `from_same_b` makes the
cc contribution small, which is a result, not grounds to drop the sample.
**Two obligations, both required:** the RDF normalizes **each mode separately to its own exclusive
cross section** (per-sample `N_gen`, never a shared denominator — that would give their N-weighted
average, ≈2× under-normalizing each), and **every plotting code must ADD the two**. A curve showing
only `bb` is not "POWHEG". Details: `docs/powheg.md`; plotting rule:
`.claude/conventions/atlas-plotting.md`.

**Powheg fullsim is obsolete** for this analysis: only a Run 2 (pp17) sample
exists, and Powheg has no pTHat slicing, so its high-p_T statistics are too poor
for high-p_T efficiency / response. It is not requested for Run 3. (Details:
`docs/powheg.md`.)

## 6. Backgrounds & purity

- **Combinatorial:** estimated from SS pairs (§4a).
- **Gluon splitting g → QQ̄:** the dominant *correlated* background — two muons
  from a single gluon-split heavy-quark pair also give low mass / small opening
  angle, mimicking the single-b signal. This must be modeled and separated from
  the single-b signal.
- **Fake muons** (π/K decay-in-flight, mis-ID): assessed via a **Δp/p
  significance template fit** (method adapted from the Run 2 note); framework to
  be built.

**Level of the template fits — RECONSTRUCTED.** Both the low-mass minv template fit and
the Δp/p fake-muon fit are performed on **reconstructed quantities**, after
trigger-efficiency correction but **before** reconstruction-efficiency correction and
the bin-by-bin detector correction (§4c); the detector corrections and signal acceptance are applied to the extracted
signal yield **after** the fit. The mixed-event combinatoric template is built from
real-data (reconstructed) muons and the fake background has no truth analogue, so the fit
cannot live at truth. (Detail: `docs/tracking/low_mass_dimuon_template_fit.md` §3.)

**Signal vs background templates (procedure under development — see
`docs/powheg.md`).** Powheg-HQ truth provides the **NLO-hard** part of the
spectrum but is **strongly suppressed in low-mass / low-k_T gluon splitting**
(kinematically, stacking down the chain: the splitting gluon needs p_T ≳ 2 m_Q to
make the QQ̄ pair; each heavy quark hadronizes and decays semileptonically to a
muon carrying only a fraction of the hadron momentum; and each of those muons
must still pass p_T > 4.5 GeV — so the gluon must be well above 2 m_Q). At NLO the classic flavour-creation / flavour-excitation /
gluon-splitting "topologies" are no longer well-separated — they merge into the
2→3 matrix element — *except* that the low-mass g→QQ̄ tail is the missing piece. The
anticipated approach is to **stitch** the Powheg-HQ NLO piece together with a
**Pythia gluon-splitting template** that supplies the low-k_T/low-mass g→QQ̄
contribution, for the background (and possibly signal-mix) templates. Exact
stitching and normalization are to be developed.

## 7. Conventions (quick reference)

- Trees: `..._sign1` = same-sign (SS), `..._sign2` = opposite-sign (OS).
- "HI" in branch names = Heavy Ion (Pb+Pb only); never apply to pp.
- Triggers: PbPb = single-mu4; pp24 = 2mu4. **mu4_mu4noL1 is not currently used**
  for either system.
- **Vertices differ between the two systems, on purpose.** Pb+Pb has ~10⁻³ collisions per
  bunch crossing AND rejects pile-up at event level, so its pairs are primary-vertex pairs:
  both muons must pass |d₀| < 2 mm and |z₀ sin θ| < 2 mm w.r.t. **the primary vertex**. pp24
  has ~4 collisions per crossing and its **luminosity is defined over all of them**, so
  restricting the numerator to the primary vertex while dividing by an all-collision L_int
  would bias dσ/dX low: pp pairs are therefore selected against **every** reconstructed
  track-bearing vertex, and are kept when both muons pass the same cuts w.r.t. **the same**
  vertex, primary or secondary (2026-09-08). Applied at the NTuple stage and **mirrored in the
  pp-conditions fullsim MC**, so ε_reco is measured on the selection it corrects.
  See `docs/tracking/pp24_all_vertex_pairs.md`.
- Reco/trigger efficiencies and all derived numbers depend on external inputs
  (lumi, ⟨T_AA⟩, σ_PbPb, GRLs, HLT names) that are **not yet confirmed**; see
  roadmap §Q2 before quoting any number.
