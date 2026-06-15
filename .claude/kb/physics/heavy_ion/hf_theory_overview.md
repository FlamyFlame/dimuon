# Open Heavy Flavor in Nuclear Collisions — Theory Overview

**Source:** P. B. Gossiaux (SUBATECH, Nantes), "Open Heavy Flavors in Nuclear Collisions: Theory Overview" — Quark Matter 2018 proceedings, *Nuclear Physics A* (2019).
**arXiv / DOI:** arXiv:1901.01606v1 [nucl-th] (6 Jan 2019)
**URL:** https://arxiv.org/abs/1901.01606
**PDF:** not committed (SUPPORTIVE — read on demand: `curl -sSL https://arxiv.org/pdf/1901.01606 -o /tmp/x.pdf`, then `gs`)
**Classification:** SUPPORTIVE / complementary (theory-modeler's companion to the experimental/big-picture HI reviews already in the KB)
**Added:** 2026-06-14

Short (8-page) **theory-modeler's** overview of heavy-quark (HQ) transport and
energy loss in the QGP. It complements the broader reviews in the KB by giving
the **modeling-side vocabulary** an agent needs when paraphrasing how our
measured **R_AA** (and a possible **v₂**) connect to medium transport
coefficients. Value is entirely **`[background-for-writing]`** — it is NOT a
source of methods or numbers we reproduce. For the shared concepts (dead cone,
collisional vs radiative loss, R_AA / v₂ definitions, spatial diffusion
coefficient D_s, the flow "bump") this doc does **not** re-explain — see
[[physics/heavy_ion/open_hf_production]]; here we keep only Gossiaux's
*theory-specific* additions.

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Fokker–Planck / Langevin transport vocabulary: drag η_D (inverse relaxation time), transverse diffusion κ_T, q̂ = 2κ_T, spatial diffusion D_s, relaxation time τ_relax ∝ m_Q (Eqs. 1–2, p.1–2) | IntNote §1 / thesis: phrasing of *what* b-quark R_AA constrains about the medium | [background-for-writing] |
| Classification of HQ-in-QGP models: transport-coefficient/Langevin vs cross-section/Boltzmann; elastic vs elastic+radiative vs radiative (Table 1, p.2) | Thesis ch.1 framing of the model landscape our R_AA is compared against | [background-for-writing] |
| The "R_AA–v₂ puzzle" and the ingredients invoked to resolve it (hadronic rescattering, low-pT coalescence, T-dependent drag; Boltzmann vs Langevin) (§2, §4) | Interpretation/caveats if we measure both b-dimuon R_AA and v₂ (analysis_overview §3, roadmap §Q2.9) | [background-for-writing] |
| Collective extraction efforts (EMMI/Beraudo brick & bulk comparisons, Duke Bayesian D_s, X.-N. Wang "brick problem") and their stated uncertainty budget (§2, §4) | Thesis context on theory uncertainty; tempering of model-to-data comparison of our R_AA | [background-for-writing] |
| Hadronization prescriptions: dual fragmentation+recombination; IPC / resonance recombination / in-medium fragmentation; ~25% model spread in H_AA (§4) | Caveat: low-/intermediate-pT b-hadron spectra (hence our pair pT) carry hadronization-model dependence | [background-for-writing] |
| Bottom-specific note: precise **v_n(B)** measurements "expected in the next runs"; mass→smaller v_n/ε_n (HQ inertia) (§3) | Motivation that a Run-3 **b** flow/suppression measurement is a wanted input | [background-for-writing] |

## Scope & condition-difference warnings
- **Theory proceedings, 2018**, mostly **charm (D-meson) R_AA / v₂** at **RHIC
  Au+Au 200 GeV and LHC Pb+Pb 2.76 / 5.02 TeV**. Bottom data are sparser; v_n(B)
  is anticipated, not yet measured here.
- Our analysis is **open-b → single-b dimuon at Run 3 5.36 TeV**, a different
  channel and energy. Treat all of this as **interpretive background only** —
  none of its numbers (D_s values, uncertainty percentages) are inputs to our
  result. ACKNOWLEDGE the energy/channel difference; DO NOT assume the size or
  direction of any consequent change (§5).
- The numerical statements below are the **models'** quoted spreads, not
  measurements and not values we adopt.

## Content summary   (theory-specific slice; relevance-first)

**Transport coefficients (Eqs. 1–2, p.1–2).** A HQ undergoing elastic scattering
is described by stochastic (Langevin) equations whose moments give a **drag**
coefficient A ≡ η_D (units fm⁻¹, interpreted as inverse relaxation time) and a
**transverse diffusion** coefficient κ_T (GeV²fm⁻¹), with the jet transport
parameter **q̂ = 2κ_T**. At low momentum (p ∼ m_Q, diffusive regime) a
generalized Einstein relation links them; for p→0, κ = 2T m_Q η_D, and one
writes everything via the **spatial diffusion coefficient D_s**:
(2πT)D_s = 4πT³/κ = 2πT²/(E_Q η_D), with τ_relax = η_D⁻¹. lQCD gives
**2πT D_s ≈ 6 ± 2 near T_c** (large uncertainty) → τ_relax ≈ (3 ± 1.5) fm.
⟨These are model/lQCD numbers, not ours.⟩

**Model landscape (Table 1, p.2).** Two families: **transport-coefficient based
(Langevin)** — e.g. TAMU, Catania-LV, POWLANG (HTL/lQCD), Duke, ASW, AdS/CFT,
DABMOD; and **cross-section / |M|²-based (Boltzmann)** — e.g. AMPT, MC@sHQ
(elastic and +rad), BAMPS, PHSD, Catania-BM, CUJET3, LIDO, Djordjevic et al.,
SCET_G,M. Columns split **elastic**, **elastic+radiative**, **radiative**, and
"other". Radiative loss adds path-length dependence Δp_rad ∝ L^α (α ≥ 2) and
coherence effects; most schemes still relate it back to the Fokker–Planck
coefficients.

**R_AA–v₂ puzzle (§2).** Qualitatively R_AA and v₂ are understood (high-pT
radiative suppression with recovery at very high pT; low-pT equilibration/flow
"bump"), but many models — especially elastic+radiative ones — predicted **too
small v₂**, and D_s of models that fit the data spanned a **factor ~5**. The
Catania group attributes the resolution to extra ingredients that build v₂ late
(near/below T_c) while R_AA is set early: (a) hadronic rescattering (~+1% v₂),
(b) low-pT coalescence (~+1%), (c) T-dependence of the drag (up to +3% going
from pQCD η_D∝T² to non-perturbative η_D∝T⁰, e.g. T-matrix/QPM/PHSD). A **full
Boltzmann transport yields ~1% larger v₂ than Langevin** at fixed R_AA (Langevin
FDT suppresses longitudinal fluctuations). Increasing evidence that **beyond-LO
pQCD** effects near T_c are needed, favoring **small D_s**.

**Coefficient-extraction efforts (§2, §4).** Duke **Bayesian** fit of D_s(T,p)
from an extended data set finds significant non-perturbative contribution up to
pT ≈ 20 GeV/c, but the T-dependence of D_s is not yet well constrained. EMMI/
Beraudo **model-comparison** (common pQCD×5 energy loss) found c-quark R_AA in
0.3–0.4 (0–10%) / 0.4–0.6 (30–50%); the **brick problem** (uniform medium, tuned
to a fixed R_AA at 15 GeV/c) shrinks the extracted-coefficient spread and
separates interaction classes (QPM / pQCD-like / elastic+radiative). Combined
theory uncertainty per model (fixed energy loss) is at least experiment-sized:
**≈15% on R_AA(D), ≈1% on v₂(D)**.

**Hadronization (§4).** Consensus on a **dual** mechanism — fragmentation at high
pT, **recombination/coalescence** at low pT — but the recombination prescription
varies: instantaneous parton coalescence (IPC, most common), resonance
recombination, in-medium fragmentation. Quantified via H_AA ≡
(dN_D/dpT)/(dN_c,final/dpT); model deviations **~25%**. Λ_c/D₀ enhancement at
intermediate pT is a recombination signature.

**High-pT & jet-framework schemes (§2).** DGLV, higher-twist, SCET_{G,M} handle
coherence and mass in radiative loss; DREENA-B adds Bjorken bulk evolution;
multi-stage (medium-modified PYTHIA + on-shell transport) reproduces D, B R_AA
for pT > 10 GeV/c. Benchmark ΔE(L) still disagrees between schemes.

**Flow & new observables (§3, §5).** v_n/ε_n decreases for heavier mass / higher
harmonic / more peripheral (HQ inertia). **v_n(B)** precision is expected in
upcoming runs. Proposed extra probes: HF azimuthal correlations, HF momentum
imbalance, directed flow v₁(D/D̄), Λ_c R_AA — mostly dominated by initial-state
effects (≲10% sensitivity to energy-loss type), so they demand high precision.

## References worth future reading   (≤3)

1. **A. Beraudo et al. (EMMI), "Extraction of Heavy-Flavor Transport Coefficients
   in QCD Matter", arXiv:1803.03824** [ref 5] — *SUPPORTIVE*. *New info:* the
   quantitative model-to-model comparison (initial spectra, bulk, transport,
   hadronization) and the H_AA hadronization-uncertainty quantification this
   overview summarizes. *Serves:* theory-uncertainty framing for our b R_AA
   interpretation/systematics discussion. (Same EMMI task force already queued
   from `open_hf_production` ref 68 — confirm not a duplicate before adding.)
2. **F. Prino, R. Rapp, "Open Heavy Flavor in QCD Matter and in Nuclear
   Collisions", J. Phys. G43 (2016) 093002, arXiv:1603.00529** [ref 4] —
   *SUPPORTIVE*. *New info:* a dedicated longer review of HQ transport including
   the D_s-vs-drag-momentum-dependence point (POWLANG-HTL vs TAMU) used here.
   *Serves:* thesis ch.1 transport-coefficient background.

*(Ref [2] Andronic et al. EPJC76 (2016) 107, arXiv:1506.03981 is already queued
via [[physics/heavy_ion/open_hf_production]] — not re-added here.)*

## Related KB docs   (knowledge graph)
- [[physics/heavy_ion/open_hf_production]] — the anchor open-HF review (Dong–Lee–Rapp); owns the shared concepts (dead cone, collisional/radiative, D_s, R_AA/v₂ definitions). This doc adds the modeling/transport-classification specifics.
- [[physics/heavy_ion/hf_hot_qcd_matter]] — HF-in-QGP review (hot-QCD-matter implications); companion experimental/phenomenology overview.
- [[physics/heavy_ion/hi_big_picture]] — Busza–Rajagopal–van der Schee QGP big-picture (R_AA, flow, Glauber motivation).
- [[analysis/overview]] — our R_AA / single-b dimuon goal; this source supplies the transport-theory background for its Intro.
