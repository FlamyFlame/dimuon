# Open Heavy-Flavor Production in Heavy-Ion Collisions

**Source:** X. Dong, Y.-J. Lee, R. Rapp, "Open Heavy-Flavor Production in Heavy-Ion Collisions" — Annu. Rev. Nucl. Part. Sci. (review)
**arXiv / DOI:** arXiv:1903.07709 [nucl-ex] (2019); 10.1146/annurev-nucl-101918-023806
**PDF:** `./1903.07709.pdf` (PRIMARY)
**Classification:** PRIMARY (key field review; the exact subfield this analysis studies)
**Added:** 2026-06-14

This is the anchor background review on **open heavy flavor (HF) in heavy-ion
collisions** — the subfield of our Run-3 single-b → dimuon measurement. Self-
sufficient for intro/motivation writing and for interpreting our R_AA; the PDF
is the fallback for figure-level model details.

---

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| HF as a "calibrated probe" of the QGP; b/c produced only in initial hard scatterings, thermalization time ~ fireball lifetime (p.3) | IntNote §1 Introduction; thesis ch.1 motivation | [background-for-writing] |
| R_AA definition + N_coll = T_AA·σ_pp^inel scaling logic; R_AA=1 = no medium effect (Eq.1, p.4) | Our primary observable R_AA (analysis_overview §3); IntNote intro framing of normalization | [background-for-writing] / interpretation |
| v_n / v_2 definition as Fourier coefficients vs event plane; two origins of v_2 (low-pT pressure-gradient flow vs high-pT path-length) (Eq.2, pp.4–5) | Potential v_2 observable (analysis_overview §3.3, roadmap §Q2.9); thesis flow background | [background-for-writing] |
| Energy-loss mechanisms: collisional (elastic) vs radiative (gluon) and where each dominates in pT (§4.1) | Interpreting our b-quark R_AA(pair pT); IntNote intro on energy loss | [background-for-writing] / interpretation |
| Dead-cone effect: gluon radiation suppressed for Θ < m_Q/E → mass hierarchy ΔE_b<ΔE_c<ΔE_q<ΔE_g (p.3, §4) | Core physics motivation for using **b** (heaviest) as probe; thesis ch.1 | [background-for-writing] |
| Bottom-charm R_AA hierarchy: non-prompt J/ψ & non-prompt D0 R_AA above prompt D0 at intermediate pT; b-c degeneracy onset ~25–30 GeV (§4.2, Fig.5) | What our single-b R_AA result should be compared to / contextualized against | [background-for-writing] / interpretation |
| HF-triggered HF–HF angular correlations as elastic-vs-radiative discriminator; distinct from HF–light v_2 (§4.3) | Context for any Δφ correlation-width measurement (roadmap §Q2.9 open question) | [background-for-writing] |
| Cold-nuclear-matter (CNM) effects in p+A: shadowing/nPDF, Cronin, final-state energy loss; small at mid-y, larger at forward y (§6.1) | Caveat when interpreting R_AA as purely hot-medium; thesis/IntNote systematics framing | [background-for-writing] |

---

## Scope & condition-difference warnings

- System/energy: review covers **RHIC Au+Au √s_NN=200 GeV** and **LHC Pb+Pb
  √s_NN=5.02 TeV** (plus p+A small systems). Our analysis is **LHC Pb+Pb 5.36 TeV
  Run 3**. ACKNOWLEDGE the energy difference (5.02→5.36 TeV); DO NOT assume the
  size/direction of any consequent change in spectra, T_AA, or suppression.
- Channel: the review centers on **D-meson** R_AA/v_2 (charm); bottom data are
  sparser (non-prompt J/ψ, non-prompt D0, B±). Our signal is **single-b → dimuon**
  (two muons from one b-hadron chain), not measured directly here — treat the b
  results as comparison/context, not as a method we reproduce.
- The review is **2019**; it pre-dates Run 3. Numbers it quotes are RHIC/Run-2 era.

---

## Content summary   (relevance-first)

### Big picture (§1)
- QGP forms above T_pc ≈ 155–160 MeV (lQCD). The medium is **strongly coupled**:
  small η/s ≈ 0.1–0.3 (near the 1/4π bound) from light-hadron flow fits; high-pT
  hadron suppression gives jet transport q̂/T³ ≈ 3–6 at ~10 GeV parton energy.
- Heavy quarks are "heavy" because (i) m_Q ≫ Λ_QCD → pQCD-calculable production,
  and (ii) m_Q ≫ T_QGP → produced **only in initial hard scatterings** (formation
  time ~1/2m_Q ~ 0.1 fm/c), thermalization time ~ fireball lifetime. So HF carries
  the medium's history across the full pT range — a "conserved" probe. >95% of HQs
  hadronize into open-HF hadrons (rest → quarkonia / multi-HF).
- pT regimes (the organizing theme): **low pT** → Brownian motion / diffusion
  (elastic, q²~T²), gauge = v_2 + low-pT R_AA "flow bump"; **intermediate pT
  ~5–10 GeV** → transition (few scatterings) + hadronization (recombination);
  **high pT >10 GeV** → energy loss, radiative-dominated.

### Observables (§1, Eqs. 1–2)
- **R_AA(pT) = [d²N_AA/dpT dy] / [N_coll · d²N_pp/dpT dy]**, with
  N_coll = T_AA · σ_pp^inel. R_AA = 1 means binary-collision scaling / no medium
  effect; energy loss pushes R_AA < 1. (Matches our analysis_overview §3 R_AA with
  ⟨T_AA⟩ normalization.)
- **v_n = ⟨cos(n(φ−Ψ_EP))⟩**, Fourier coefficients of azimuthal distribution vs
  event plane. v_2 ("elliptic flow") is largest. Two sources: (a) low-pT pressure
  gradients in the thermalizing almond-shaped overlap → positive v_2 (probes
  viscosity); (b) high-pT path-length difference (shorter in-plane path → less
  suppression) → positive v_2 (probes energy-loss path dependence).

### HQ diffusion at low pT (§2)
- Framework: Fokker-Planck / Langevin; the key transport coefficient is the
  **spatial diffusion coefficient D_s** (D_s = T / [m_Q A(p=0)]), analogous to
  η/s and electric conductivity. The dimensionless **2πT D_s** carries universal
  QGP transport info.
- D-meson data (RHIC + LHC) show large v_2, peaking **>15% at pT≈3 GeV/c**, falling
  to ~5% above 10 GeV/c, with a corresponding low-pT R_AA "flow bump" → clear
  charm collectivity. LHC and RHIC R_AA agree within uncertainties despite ~25×
  energy difference.
- Phenomenology converges on **2πT D_s ≈ 2–4 near T_pc** (charm); pure LO-pQCD
  values are ruled out. Most models prefer an increasing D_s with temperature.
  Radiative-included models (Nantes/MC@sHQ, BAMPS) need smaller D_s than
  elastic-only (TAMU, PHSD) because 3-body radiative final states generate less v_2.

### Hadronization (§3)
- High pT: fragmentation functions (universal). Low/intermediate pT: **coalescence
  / recombination** with comoving QGP partons; signatures = baryon/meson
  enhancement, constituent-quark-number scaling (CQNS) of v_2.
- HF hadro-chemistry probes hadronization: **D_s+/D0** enhanced in A+A (~0.35,
  consistent with statistical hadronization), **Λ_c+/D0** enhanced by ≥2× over pp,
  hint of enhanced **B_s0/B+** (CMS). Recombination clearly at work. Large Λ_c and
  D_s contributions mean total charm cross section ~ N_bin scaling once all states
  summed.

### Energy loss at high pT (§4) — most relevant to our b probe
- "Jet quenching": hard partons lose energy traversing the QGP, suppressing high-pT
  R_AA (first at RHIC, confirmed LHC). Two mechanisms: **collisional** (elastic
  scattering off medium) and **radiative** (medium-induced gluon emission;
  radiated gluon can re-scatter).
- **Dead cone:** gluon radiation off a heavy quark is suppressed at angles
  Θ < m_Q/E → less radiative loss for heavier/lower-E quarks → expected hierarchy
  **ΔE_b < ΔE_c < ΔE_q < ΔE_g** → less HF suppression than light hadrons. (Gluons
  lose most due to larger color charge.)
- **D-meson R_AA shape (Fig.2):** drops to a minimum at **5–8 GeV/c**, then rises,
  reaching ~0.7 at 100 GeV/c (CMS). Low-pT region = diffusion; minimum ≈ onset of
  radiative loss; high-pT rise driven by flattening production spectrum + reduced
  coupling.
- **Mass hierarchy (Fig.5, §4.2):** at intermediate pT, non-prompt J/ψ R_AA (from
  b→J/ψ) and non-prompt D0 R_AA (b→D0) sit **above** prompt-D0 R_AA → bottom less
  suppressed than charm. (Caveat: decay products carry less pT than parent; a
  +~2 GeV/c horizontal shift to the b-decay R_AA strengthens the difference.) At
  very high pT, mass effects cease: CMS shows **charm-light degeneracy for
  pT>10 GeV/c**, and a further **bottom degeneracy for pT≳25 GeV/c** — consistent
  with the ~3× b/c mass ratio at the same Lorentz-γ. RHIC b→e vs c→e electrons show
  a similar hierarchy with larger uncertainties.
- Quantitative flavor conclusions need data-model comparison: primordial pT slopes
  and elastic/radiative partition differ between b and c. Path-length dependence
  ≈ linear for bottom up to ~40 GeV/c; light flavors steepen first, merging with
  c,b near 100 GeV/c at ~L^1.4 (§5.3).
- **§4.3 HF–HF angular correlations** proposed as an extra elastic-vs-radiative
  discriminator (distinct from HF–light v_2): tests few-large-Q vs many-small-Q
  deflections. So far mostly p+p (D-D, D0-hadron); first D0-jet correlation in
  Pb+Pb (CMS) consistent with HQ diffusion.

### Small systems / CNM (§6)
- p/d+A used to isolate **cold-nuclear-matter** effects: (a) nuclear-modified PDFs
  (shadowing) affecting initial HF production; (b) Cronin broadening / initial-state
  energy loss; (c) final-state interactions in the nuclear remnant. These factorize
  poorly but must be separated from hot-medium effects in A+A.
- HF R_pA at mid-y is **small** (slightly suppressed at low pT ~ shadowing,
  compatible with 1 above); larger forward-y suppression (low Bjorken-x). Cleaner
  nPDF constraints come from EW bosons, dijets, UPC charmonium.
- Surprise: **large v_2 of J/ψ, D, HF leptons in high-multiplicity p+A** (almost as
  large as light flavor) while R_pA ≈ 1 — seemingly conflicting; origin of HF flow
  in small systems not understood.

### Conclusions (§7) — summary points
1. HF is a conserved probe spanning low-pT diffusion to high-pT energy loss.
2. D-meson R_AA + v_2 → strong charm coupling; **2πT D_s ≈ 2–4 near T_pc**
   (2–4× quantum bound), scattering rates ≳ 2–3/(fm/c).
3. D-meson data show a **diffusion→energy-loss transition over pT≈5–10 GeV/c**;
   above it, D-light degeneracy; **bottom merges into degeneracy at pT≈20–30 GeV/c**,
   consistent with b/c mass ratio.
4. HF hadro-chemistry (D_s+/D0, B_s0/B+, Λ_c+/D0 enhancements) → recombination with
   thermal QGP partons.
5. Small systems: large v_2 but R_AA≈1 → conflicting QGP-formation indications.

---

## References worth future reading   (≤3)

1. **A. Andronic et al., "Heavy-flavour and quarkonium production in the LHC era…",
   Eur. Phys. J. C76:107 (2016)** [ref 37] — *SUPPORTIVE*. New info: a far more
   detailed, dedicated review of HF + quarkonium production from pp to A+A
   (pp baseline, production mechanisms, energy-loss models). Serves IntNote §1
   Introduction and thesis ch.1 background more thoroughly than this survey.
2. **A. Beraudo et al. (EMMI task force), Nucl. Phys. A979:21 (2018)** [ref 68] —
   *SUPPORTIVE*. New info: quantitative model-to-model comparison of HF energy loss,
   path-length dependence, and the elastic/radiative partition (incl. bottom).
   Serves interpretation of our b-quark R_AA and the energy-loss systematics
   discussion.
3. **CMS bottom R_AA measurements (B± and non-prompt J/ψ / non-prompt D0),
   Sirunyan et al. — behind Fig. 5** [refs 118, 120] — *PRIMARY (comparison data)*.
   New info: actual measured **bottom** R_AA(pT) at LHC, the closest existing
   measurements to compare our single-b R_AA against. Serves IntNote Results
   comparison and thesis prior-measurement context.

---

## Related KB docs   (knowledge graph)

> **This is the canonical HF-physics hub** (energy loss, R_AA, dead-cone, transport
> coefficients). The other heavy-ion docs defer their shared concepts here.

- [[hi_big_picture]] — broader QGP/HIC frame (η/s, jet quenching, Glauber) this HF review sits inside.
- [[hf_hot_qcd_matter]] — earlier (2013) open-HF review; RHIC-era discovery + decay-lepton/DD̄-dimuon detail.
- [[rhic_open_hf_review]] — RHIC-energy experimental complement (open-b via DCA, RHIC↔LHC contrast).
- [[hf_theory_overview]] — theory-modeling companion (Langevin/Boltzmann, transport-coefficient extraction).
- [[../background/gluon_splitting_flavour_excitation]] — the g→QQ̄ / DD̄ correlated-background mechanism.
- [[analysis/overview]] / `Analysis/docs/analysis_overview.md` — our R_AA / single-b dimuon goal (this doc is its field background).
- [[analysis/run2_dimuon_note]], [[analysis/run2_hf_muon_raa]], [[analysis/run2_hf_muon_vn]] — the prior measurements this physics motivates.
</content>
</invoke>
