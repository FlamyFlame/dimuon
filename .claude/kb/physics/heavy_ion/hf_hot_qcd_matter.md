# Heavy-Flavor Production in Heavy-Ion Collisions (Averbeck review, 2013)

**Source:** R. Averbeck, "Heavy-flavor production in heavy-ion collisions and implications for the properties of hot QCD matter" — Prog. Part. Nucl. Phys. **70** (2013) 159–209 (review; experimental results published until July 2012).
**arXiv / DOI:** doi:10.1016/j.ppnp.2013.01.001 (no arXiv id on title page).
**PDF:** `./Heavy-flavor production in heavy-ion collisions and implications for the properties of hot QCD matter.pdf` (PRIMARY)
**Classification:** PRIMARY (foundational, heavily-citable PPNP field review of open HF in HIC).
**Added:** 2026-06-14

The **earlier/foundational** review of the exact subfield of our analysis (open
heavy flavor in heavy-ion collisions), covering the **RHIC-era establishment**
plus the first LHC Pb+Pb results. Value is almost entirely
**`[background-for-writing]`** for the IntNote §1 Intro and thesis ch.1.

> **DEDUP — read this first.** The modern, more relevant review of this subfield is
> already in the KB: **[[heavy_ion/open_hf_production]]** (Dong/Lee/Rapp 2019,
> LHC 5.02 TeV era). That doc owns the **modern transport-coefficient** content —
> spatial diffusion coefficient **2πT D_s ≈ 2–4 near T_pc**, the diffusion→energy-loss
> pT transition, the bottom/charm degeneracy onset (~20–30 GeV). **This 2012-era
> review does NOT use the 2πT D_s formulation** (it speaks of heavy-quark diffusion
> giving "indirect access to η/s" only). Do not quote a 2πT D_s number from here.
> Use this Averbeck doc for what it adds *beyond* Dong/Lee/Rapp: the RHIC-era
> discovery narrative, the **semileptonic (e/μ) decay-lepton channel** detail, and
> the **intermediate-mass-region DD̄ dimuon continuum**.

---

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| HF (c,b) produced **only** in initial hard pQCD scatterings → total yield scales with N_coll (binary scaling); deviations = nuclear modification; thermal HF production negligible even at LHC (§1.2, p.163–164) | IntNote §1 Intro; justifies HF as a calibrated probe | [background-for-writing] |
| Energy-loss mass hierarchy **R_AA^π < R_AA^c < R_AA^b** from dead-cone suppression of radiative loss; collisional (elastic) loss comparable and *more* important for heavy quarks (§1.2, p.164) | Core motivation for using **b** as probe; interpreting our b R_AA(pair pT) | [background-for-writing] / interpretation |
| Heavy-quark transport models (Langevin diffusion → indirect η/s; AdS/CFT drag; in-medium D/B resonance dissociation) confronted with data (§1.2 p.164; §4.3 p.202) | Intro framing of "what HF teaches about QGP transport"; defer numbers to [[heavy_ion/open_hf_production]] | [background-for-writing] |
| **v_2 of HF decay leptons** sensitive to heavy-quark thermalization; non-zero charm v_2 at RHIC → charm flows; D0 v_2 hint at LHC (§1.2 p.165; §3.4 p.191; §4.3 p.202) | Motivation for the **potential v_2 observable** (analysis_overview §3.3, roadmap §Q2.9) | [background-for-writing] |
| **Semileptonic (e/μ) decay-lepton method**: measure inclusive e/μ, subtract cocktail/decay-muon background, R_AA & v_2 of HF leptons (§3.4, §4.1.1) | Closest measured analog to our **muon-based** measurement; background-subtraction framing | [background-for-writing] |
| **ALICE forward-muon R_AA** (c,b→μ, 2.5<y<4.0, Pb+Pb 2.76 TeV): strong, flat suppression in central, weaker peripheral; high-pT muons have **significant bottom fraction** → bottom also strongly coupled (§4.3.2 p.202–203) | The single most directly comparable prior measurement to our single-b→dimuon R_AA | [background-for-writing] / comparison |
| **Intermediate-mass-region (IMR) dimuon continuum** from correlated **DD̄** semimuonic decays; in-medium energy loss + opening-angle modification; at LHC bottom contributes significantly (§2.1 p.166; §5 p.204) | Background framing: our **g→QQ̄** background = two muons from *two* heavy quarks; any Δφ/opening-angle discussion | [background-for-writing] |
| FONLL / GM-VFNS NLO pQCD describe HF pp production; PYTHIA needs K-factor 2–4.5 (charm), ~2 (bottom) (§1.2 p.163; §2.1 p.165–166) | pp-baseline framing; context for our Powheg/Pythia template choices (analysis_overview §6) | [background-for-writing] |

---

## Scope & condition-difference warnings

- Experimental data reviewed are **RHIC Au+Au/Cu+Cu √s_NN = 0.2 TeV** and the
  **first LHC Pb+Pb √s_NN = 2.76 TeV** results (published ≤ July 2012). Our
  analysis is **Run 3 Pb+Pb 5.36 TeV**. ACKNOWLEDGE the large energy/era gap; DO
  NOT port any quoted suppression magnitude, energy density, or T.
- Channel: results are HF **decay electrons/muons** (mid- and forward-rapidity)
  and **D mesons**; bottom is mostly inferred indirectly. Our signal is
  **single-b → dimuon** (two muons from *one* b-hadron chain) — not measured here;
  treat all results as comparison/context, not a method we reproduce.
- This is a **2013** review — pre-dates the modern transport-coefficient
  consensus; for current numbers use [[heavy_ion/open_hf_production]].

---

## Content summary   (relevance-first; distinctive content only)

### HF as a QGP probe (§1.2)
Charm/bottom masses (m_c≈1.29, m_b≈4.19 GeV) ≫ Λ_QCD → pair production is a hard
pQCD process at *all* momenta (the mass is a hard scale even at p=0). Consequences
used in our Intro: (1) HF pp production tests pQCD and provides the AA baseline;
(2) total HF yield is set by initial hard scattering → **binary (N_coll) scaling**
in AA, broken only by initial-state nPDF effects, not by final-state energy loss
(which redistributes pT but conserves yield); (3) thermal HF production in the
fireball is negligible. **Initial-state effects** (shadowing at low x, anti-shadowing
at x~0.1, gluon saturation) must be disentangled from **final-state** (hot-medium)
effects.

### Energy loss & the mass hierarchy (§1.2 p.164) — core motivation
- **Radiative** loss (medium-induced gluon radiation) is the historically dominant
  mechanism; it is **gluon > light-quark** (color factor) and **suppressed for
  heavy quarks** by the **dead-cone effect** (small-angle radiation cut off by the
  quark mass) → expected ordering at pT≲10 GeV: **R_AA^π < R_AA^c < R_AA^b**.
- **Collisional** (elastic) loss is comparable in magnitude over a wide kinematic
  range and is **relatively more important for heavy quarks** — RHIC data forced
  this realization (radiative-only models could not explain the large HF-electron
  suppression).
- Other Ansätze: Langevin transport (diffusion + damping → indirect η/s), repeated
  in-medium HF-hadron formation/dissociation, AdS/CFT string-drag. The **ratio
  R_AA^c/R_AA^b is predicted to differ markedly between AdS/CFT and pQCD** (Fig.27)
  → motivates *separate* charm/bottom measurements. (Modern numeric outcome: see
  [[heavy_ion/open_hf_production]].)

### R_AA & v_2 as the observables (§1.2)
- Departures of HF spectra from N_coll-scaled pp = the energy-loss signature
  (softening of the pT spectrum). For light hadrons R_AA reaches its minimum at
  pT≈6–7 GeV/c and stays <1 up to 100 GeV/c at the LHC.
- **v_2** (elliptic flow, Fourier coefficient vs reaction plane) probes the degree
  of heavy-quark **thermalization**/collective participation; non-zero HF v_2 ⇒
  strong coupling. Simultaneous R_AA + v_2 is the more stringent model test.

### RHIC-era establishment (§3.4) — the discovery narrative
- HF **decay-electron** R_AA at mid-rapidity in central Au+Au is suppressed
  ~as strongly as light hadrons — surprising vs radiative-only expectations →
  collisional loss / alternative mechanisms required. Measured electron **v_2 ≠ 0**
  ⇒ charm participates in collective flow. Together these branded the QGP a
  **strongly coupled, near-perfect fluid** (minimal η/s). Caveat stated by the
  author: semileptonic inclusive electrons give only **indirect** info — charm/bottom
  cannot be cleanly separated without vertex detectors.
- HF **decay-muon** R_AA at forward rapidity (PHENIX Cu+Cu) shows suppression
  comparable to mid-rapidity Au+Au electrons — surprising given lower forward
  energy density → possible cold-nuclear-matter role at forward y.

### First LHC Pb+Pb (§4.3) — directly relevant channel
- **D-meson R_AA** (ALICE, 2.76 TeV, central) strongly suppressed, similar to
  charged hadrons; models combining radiative + collisional/dissociation describe
  it reasonably (tendency to under-predict). AdS/CFT-drag under-predicts.
- **D0 v_2** (ALICE) shows an indication of non-zero charm flow at the LHC (large
  uncertainties).
- **HF decay-muon R_AA** (ALICE, 2.5<y<4.0): strong, **pT-flat** suppression in
  0–10% central over 4<pT<10 GeV/c, weaker in 40–80%; pronounced centrality
  dependence. **A significant fraction of these high-pT muons comes from bottom**
  → bottom quarks are *also* strongly coupled to the medium. This is the closest
  existing measured analog to our single-b→dimuon R_AA.

### Decay-lepton method & the DD̄ dimuon continuum (§4.1.1, §2.1, §5)
- Method (ALICE/PHENIX): reconstruct inclusive e or μ; subtract a **cocktail**
  (π0/η Dalitz, conversions, light VM, quarkonia, Drell-Yan) for electrons, or the
  **decay-muon** background (π/K decays, extrapolated from measured pion/kaon
  spectra in Pb+Pb where MC fails) for muons; the remainder = HF decay leptons.
- **IMR dimuons:** correlated semimuonic decays of a **DD̄ pair** populate the
  intermediate-mass (φ–J/ψ) dimuon continuum, alongside Drell-Yan and thermal
  radiation. In-medium effects modify both the energy loss *and* the **opening
  angle** between the two heavy quarks (full thermalization ⇒ angular correlation
  washed out). At the LHC, **bottom** contributes significantly to the IMR.
  *(Note for our work: this DD̄ two-heavy-quark dimuon source corresponds to our
  g→QQ̄ correlated background — distinct from our single-b signal where both muons
  come from one b chain; analysis_overview §6.)*

---

## References worth future reading   (§6; ≤3)

> Most high-value HF refs (Andronic 2016, Beraudo 2018 energy-loss task force, CMS
> bottom R_AA) are already queued under [[heavy_ion/open_hf_production]]; not
> repeated here. Only refs **distinctive to this review** are listed.

1. **ALICE Collaboration, HF decay-muon R_AA at forward rapidity in Pb+Pb 2.76 TeV**
   (ref [222]; ©2012 APS) — *PRIMARY (comparison data)*. New info: the actual
   measured **muon-channel** HF R_AA (the same observable/decay channel as our
   dimuon-from-b measurement), with explicit bottom contribution at high pT.
   Serves IntNote Results comparison + thesis prior-measurement context.
2. **A. Majumder, M. van Leeuwen, "The theory and phenomenology of perturbative
   QCD based jet quenching", Prog. Part. Nucl. Phys. 66 (2011) 41** (ref [90]) —
   *SUPPORTIVE*. New info: the dedicated theory/phenomenology review of parton
   energy loss (radiative vs collisional formalism, q̂) underlying the R_AA
   interpretation; deeper than any current KB doc. Serves the energy-loss
   discussion in the Intro and systematics framing.

(No further references worth adding — the DD̄/IMR-dimuon refs are SPS-era and
low-relevance, and the modern transport refs are already covered elsewhere.)

---

## Related KB docs   (knowledge graph)
- [[heavy_ion/open_hf_production]] — **the modern (2019) review of this same
  subfield**; owns the transport-coefficient (2πT D_s) numbers, the
  diffusion→energy-loss transition, and the b/c degeneracy onset. *This Averbeck
  doc is the earlier RHIC-era complement adding the decay-lepton/DD̄-dimuon detail.*
- [[heavy_ion/hi_big_picture]] — general QGP/HIC big picture (η/s, R_AA, Glauber,
  flow); the broader frame this HF review sits inside.
- [[heavy_ion/rhic_open_hf_review]] — RHIC-energy experimental complement (open-b
  via DCA template fit; RHIC↔LHC contrast).
- [[heavy_ion/hf_theory_overview]] — HF-transport theory companion.
- [[analysis/overview]] — our R_AA / single-b dimuon goal and g→QQ̄ background that
  this field background motivates.
- [[analysis/run2_dimuon_note]] — the prior Run-2 single-b dimuon measurement whose
  Intro/motivation this review supports.
</content>
</invoke>
