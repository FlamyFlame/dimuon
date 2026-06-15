# Heavy Ion Collisions: The Big Picture and Big Questions

**Source:** W. Busza, K. Rajagopal, W. van der Schee, "Heavy Ion Collisions: The Big Picture, and the Big Questions" — Ann. Rev. Nucl. Part. Sci. 2018, **68:1–49**.
**arXiv / DOI:** arXiv:1802.04801v2 [hep-ph] (2018) — doi:10.1146/annurev-nucl-101917-020852
**PDF:** `./1802.04801.pdf` (PRIMARY)
**Classification:** PRIMARY (the field's most-cited accessible QGP/HIC review; primary value is motivation text)
**Added:** 2026-06-14

Pedagogical, citation-limited review of the quark–gluon plasma (QGP) and
ultrarelativistic heavy-ion collisions (HIC). Its value to us is almost entirely
**`[background-for-writing]`** — framing *why* we study b-quark energy loss/flow
in the QGP. It is **not** a method or number source we reproduce.

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| QGP = strongly coupled near-perfect liquid, η/s ≈ 1/4π; "smallest hottest droplet of liquid made on earth" | IntNote §1 Intro / thesis ch. 1 — opening motivation | [background-for-writing] |
| Why study HIC (3 motivations: QCD cosmology / QCD phase diagram / emergence of complex matter) | Intro / thesis motivation framing | [background-for-writing] |
| **Parton energy loss & jet quenching** in QGP; dE/dx; heavy quarks as colored "hard probes" that X-ray the medium | Intro motivation for **b-quark in-medium energy loss** (analysis_overview §1) | [background-for-writing] |
| **Nuclear modification factor R_AA** definition (Eq. 5) and its interpretation as suppression | Motivation/definition for our **primary observable R_AA** (analysis_overview §3) | [background-for-writing] |
| **Glauber model**: N_part, N_coll, spectators, impact parameter, centrality classes from multiplicity percentiles | Centrality & normalization framing (analysis_overview §4d–e, roadmap §Q2.2–2.3) | [background-for-writing] |
| Charm/bottom do **not** reach chemical equilibrium → heavy-quark yields retain memory of initial hard production | Justifies heavy flavor as a calibrated probe of the medium (Intro) | [background-for-writing] |
| Azimuthal anisotropy / flow (v_n): spatial→momentum anisotropy via low-η/s hydro | Motivation for the **potential v₂ / flow** observable (analysis_overview §3, roadmap §Q2.9) | [background-for-writing] |
| Quarkonia (Υ sequential melting via Debye screening; J/ψ regeneration) | Contrast for the Intro: our signal is **open** b (single-b dimuon), not bound bb̄ | [background-for-writing] |

## Scope & condition-difference warnings
Theory/pedagogy review; experimental facts quoted are RHIC + LHC up to PbPb
**√s_NN = 2.76 / 5.02 TeV** (Run 1/2). Our analysis is **Run 3 5.36 TeV**. This
source is a motivation/conceptual reference, **not** a source of inputs or
numbers we reproduce — none of its quoted numbers (energy density, T, lumi) enter
our results. Quantities like ⟨T_AA⟩, σ_PbPb at 5.36 TeV come from elsewhere
(roadmap §Q2), not here. ACKNOWLEDGE the energy/run difference; do not port any
value.

## Content summary

**QGP, the big picture (Intro, §1; p.3–5).** Two Lorentz-contracted nuclei
collide; most interactions are soft. ~1 fm/c after impact a hot region of energy
density **far above 500 MeV/fm³** (the energy inside a hadron) forms — e.g.
**>12 GeV/fm³** at 2.76 TeV (≈20× a hadron), inferred from total transverse
energy **1.65±0.1 TeV** in |η|<0.5 (p.3). Up to **~30,000 particles** in central
PbPb (p.3). This matter is **not** hadrons and **not** a weakly coupled gas: it is
a **strongly coupled liquid** (QGP) that flows hydrodynamically with specific
viscosity **η/s ≈ 1/4π** — the smallest of any known fluid (p.3, refs 5,6; the
1/4π value is the holographic / AdS-CFT strong-coupling result, refs 30,20).
Hard scatterings at the earliest times make rare high-p_T partons → jets.

**Why study HIC (§2; p.6–13).** Three motivations: (1) **QCD in cosmology** —
HIC recreate the matter that filled the universe ~microseconds after the Big Bang
(little bangs); the QCD transition is a **continuous crossover**, not first-order
(lattice, refs 27,28). (2) **QCD phase diagram** in (T, baryon chemical potential
μ_B): crossover near μ_B=0; a possible **critical point** at higher μ_B probed by
the RHIC Beam Energy Scan; color-superconducting quark matter at high density/low
T (neutron-star cores). (3) **Emergence of complex matter** — QGP is the simplest
complex matter, closest to its QCD first principles; jets, heavy quarks of varying
mass, and quarkonia of varying size probe it across length scales.

**Lattice thermodynamics (Fig. 2; p.7).** Pressure, energy density, entropy
density show a **crossover around T ≈ 150 MeV** from a hadron resonance gas to
QGP; at accessible T they sit ~20% below the Stefan–Boltzmann (free-parton) limit
(strong coupling). The rise in ε/T⁴, s/T³ reflects color deconfinement (more
d.o.f.).

**Phenomenology & centrality (§3; p.13–21).**
- **Glauber model** (refs 53,57): nucleons that collide ≥once are **participants
  (N_part)**; those that miss are **spectators (N_spec)**; **N_coll** = number of
  binary nucleon–nucleon collisions (assuming transparency). Impact parameter *b*
  is not measurable; **centrality classes** are defined by percentiles of measured
  multiplicity/energy, mapped (monotonically) to N_part/N_coll from a Glauber MC
  (e.g. top 5% multiplicity ↔ 5% most central). Our analysis instead uses FCal ΣE_T
  → 2023 Glauber thresholds (analysis_overview §4e).
- **Hard probes scale with N_coll.** Colorless probes (γ, Z⁰) have **R_AA = 1** in
  AA — an independent validation of N_coll / the Glauber procedure (p.15). Colored
  hard probes (jets, high-p_T hadrons, **heavy quarks**) are modified by the
  medium → carry information about it.
- **Baryon stopping:** each participant loses ~2 units of rapidity (≈85% of its
  energy) in AA (p.16). **Participant scaling** of multiplicity. dN/dy ∝ s_NN^0.15.

**Nuclear modification factor (§7, Eq. 5; p.36) — our primary observable.**
```
                  dN_AA / dp_T
 R_AA(p_T) = ───────────────────────
              ⟨N_coll⟩ · dN_pp / dp_T
```
R_AA < 1 signals suppression (energy loss); R_AA = 1 for unmodified (colorless)
probes. The steeply falling jet spectrum (∝ p_T^−6) makes R_AA a sensitive
energy-loss probe (p.36). *Our analysis uses the equivalent ⟨T_AA⟩ normalization*
(R_AA = (1/⟨T_AA⟩)·(dN_PbPb/dX)/(dσ_pp/dX), analysis_overview §3) — the standard
T_AA-form of the same nuclear-modification ratio; the paper writes the N_coll form.

**Jets & energy loss in QGP (§7; p.35–39).** A parton traversing QGP (i) loses
energy/forward momentum (well established, refs 162,163), (ii) gains transverse
momentum (broadening, not yet seen), (iii) deposits energy into a **wake**. Jets
lose often **>10 GeV** over a few fm → enormous dE/dx → direct evidence the medium
is strongly coupled. dE/dx is parametrized differently for weakly- vs
strongly-coupled (holographic) plasma; data constrain its magnitude and its T-,
x-, E-dependence. Heavy quarks and high-p_T hadrons "come from jets" and are part
of this hard-probe program (ref 69). **R_AA for hadrons rises at the highest p_T**
(selection of unusually narrow jets) while jet R_AA stays suppressed — a caveat
when comparing hadron- vs jet-level suppression.

**Heavy flavor specifics (p.19–21).** Charm and bottom quark densities **do not
reach chemical equilibrium** (temperature too low) — their multiplicities
**retain memory of their initial hard production**, which is what makes heavy
flavor a calibrated probe. Quarkonia: **Υ states melt sequentially** by Debye
screening (1S survives, 2S/3S suppressed; Fig. 5, ref 99); **J/ψ regenerates** at
the LHC because ~30 cc̄ pairs/collision allow recombination of independently
diffusing c, c̄ — direct evidence quarks are deconfined. (Our signal is **open**
b → single-b dimuon, distinct from these bound-state probes.)

**Hydrodynamics & flow (§3–4; p.18–22).** Initial geometric/lumpy spatial
anisotropy of the overlap region → pressure gradients → **azimuthal momentum
anisotropies (v_n)** that survive because η/s is tiny. Measured azimuthal
correlations are well described by relativistic hydro; this is the central
evidence that QGP is a low-viscosity liquid. High-p_T v₂ exists because jets lose
less energy along the short medium axis (path-length dependence). Small systems
(pp, pPb, dAu) show similar anisotropies — open question whether proton-sized QGP
droplets form.

**Big questions (§8; p.39–41).** How/how fast QGP hydrodynamizes within ~1 fm/c;
limits of hydrodynamics / smallest droplet; how a strongly coupled liquid emerges
from an asymptotically free theory (jets as microscopes); EIC connection to
initial gluon distribution; experimental T determination / number of d.o.f.;
critical point at nonzero μ_B; first-principles multiplicity & hadronization;
color superconductivity in neutron stars.

## References worth future reading   (§6; ≤3 — open-HF weighted)

- **Connors, Nattrass, Reed, Salur**, "Jet measurements in heavy ion physics",
  Rev. Mod. Phys. **90**, 025005 (2018), arXiv:**1705.01974** (ref 69) —
  **SUPPORTIVE**. *New info:* a comprehensive **experimental** review of hard
  probes — jets, high-p_T hadrons, and **heavy-flavor (D/B, heavy-flavor-tagged)
  R_AA** measurements and methods — i.e. the actual energy-loss data this big-
  picture review only summarizes. *Serves:* IntNote §1 Intro and the R_AA
  discussion; the prior open-HF/jet-suppression measurements we will cite and
  compare against.

*Note on open heavy flavor:* this review cites **no dedicated open-heavy-flavor
review** (e.g. heavy-quark transport/diffusion coefficient, D-meson or
b-tagged-jet R_AA). That is a genuine gap to fill from elsewhere when writing the
b-energy-loss motivation — ref 69 above is the closest entry it offers.

## Related KB docs   (knowledge graph)
- [[analysis/overview]] — our physics goal, signal, R_AA observable (this source supplies the field background for that overview's Intro)
- [[analysis/run2_dimuon_note]] — the Run 2 dimuon measurement whose Intro/motivation this background supports
</content>
</invoke>
