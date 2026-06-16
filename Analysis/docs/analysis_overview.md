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
2023 / 2024 / 2025** (√s_NN = 5.36 TeV). Run 2 code paths (pp17, PbPb15/18)
exist only for cross-checks and are not maintained for new analysis decisions.

## 2. The signal: a single-b dimuon

The signal muon pair comes from **one b-hadron's decay chain**: one muon from the
direct semileptonic decay (b → μ X) and the other from the sequential charm decay
in the same chain (b → c → μ X). Because both muons share a single b parent, the
pair has characteristically **low invariant mass** and **small opening angle**.

**Truth-level signal region** (as implemented):
- opposite-sign (OS) muon pair,
- m_μμ ∈ **[1.08, 2.9] GeV**,
- pair p_T > **8 GeV**,
- |pair η| < **2.2**,
- ΔR(μ,μ) > **0.05**.

This single-b topology is what distinguishes the signal from the dominant
heavy-flavour background where the two muons come from **two different** heavy
quarks (see §6).

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
  - pp24 2mu4: ε_trig^pair = ε₁^nc · ε₂^nc.
  A pair-level ΔR trigger-correlation correction ε_dR is foreseen but currently a
  dummy ≡ 1 (measuring it needs unbiased-trigger MC).
- **ε_reco^pair** — *reconstruction* efficiency, a single **pair** efficiency
  binned in **(pair p_T, pair η, ΔR)**, measured in **Pythia fullsim** (pp) and
  **Pythia fullsim HIJING overlay** (PbPb, per centrality). It is **not** the
  product of two single-muon efficiencies: the two close-by muons are correlated
  at reconstruction (shared ID hits, ambiguity resolution), and the ΔR axis
  captures that correlation. In central Pb+Pb the efficiency is genuinely lower
  and centrality-dependent (real occupancy effect from the underlying event), not
  a bookkeeping artifact.

**(c) Detector response / unfolding.** Momentum response and bin migration are
characterized from **Pythia fullsim** and used to unfold the measured spectra to
particle level.

**(d) Normalization.** pp: integrated luminosity. Pb+Pb: number of events and the
Glauber ⟨T_AA⟩ per centrality (and the total hadronic cross-section for centrality
calibration). These external numbers must be confirmed before entering results
(roadmap §Q2).

**(e) Centrality (Pb+Pb).** From forward calorimeter (FCal) ΣE_T with 2023 Glauber
thresholds; cross-year FCal scaling brings 2024/2025 onto the 2023 scale.

**(f) R_AA & systematics.** Form R_AA from the corrected, normalized pp and PbPb
spectra; propagate efficiency, unfolding, normalization, and selection
systematics.

> The code realizes (a)–(f) through the three-stage pipeline (NTuple processing →
> RDF histogram filling → plotting). That staging is an implementation choice and
> is documented in `README.md`; it is not part of the methodology above.

## 5. MC samples and their roles

| Sample | Role |
|--------|------|
| **Pythia truth** (5.36 / 5.02 TeV) | Generator-level signal/flavour & origin studies |
| **Pythia fullsim pp24** | **Reconstruction efficiency** and **detector response / unfolding** (pp) |
| **Pythia fullsim HIJING overlay** | **Reconstruction efficiency** (PbPb), per centrality; overlay puts HIJING into the detector at hit level |
| **Powheg truth (bb / cc)** | **NLO signal template** for the template fit (see §6) |

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

**Signal vs background templates (procedure under development — see
`docs/powheg.md`).** Powheg-HQ truth provides the **NLO-hard** part of the
spectrum but is **strongly suppressed in low-mass / low-k_T gluon splitting**
(kinematically, stacking down the chain: the splitting gluon needs p_T ≳ 2 m_Q to
make the QQ̄ pair; each heavy quark hadronizes and decays semileptonically to a
muon carrying only a fraction of the hadron momentum; and each of those muons
must still pass p_T > 4 GeV — so the gluon must be well above 2 m_Q). At NLO the classic flavour-creation / flavour-excitation /
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
- Reco/trigger efficiencies and all derived numbers depend on external inputs
  (lumi, ⟨T_AA⟩, σ_PbPb, GRLs, HLT names) that are **not yet confirmed**; see
  roadmap §Q2 before quoting any number.
