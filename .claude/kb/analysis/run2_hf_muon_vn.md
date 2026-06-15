# Run 2 HF-Muon Azimuthal Anisotropy (v_n) — Pb+Pb 5.02 TeV

**Source:** ATLAS Collaboration, *Measurement of azimuthal anisotropy of muons
from charm and bottom hadrons in Pb+Pb collisions at √s_NN = 5.02 TeV with the
ATLAS detector* — Phys. Lett. B 807 (2020) 135595; CERN-EP-2019-274.
**arXiv / DOI:** 2003.03565 [nucl-ex] / 10.1016/j.physletb.2020.135595
**PDF:** `./2003.03565.pdf` (PRIMARY, committed)
**Companion internal note:** `./Run2 HF muons azimuthal anisotropies internal notes.pdf`
(ATLAS HION-2019-11, draft v0.6, Lim/Hu/Perepelitsa/Nagle/Hill, U. Colorado) —
read for implementation detail (event-plane construction, efficiency maps,
template fit). Page anchors below tagged "[PLB p.N]" or "[note p.N]".
**Classification:** PRIMARY
**Added:** 2026-06-14

---

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| **Event-plane (EP) flow-extraction method** for HF muons: FCal q-vector, recentering/flattening, two-subevent resolution Res{nΨ_n}, Fourier fit of yield vs n(φ−Ψ_n), v_n = v_n^raw/Res{nΨ_n} | The flow-fit machinery we'd need **if** v_2/v_3 of single-b dimuons is added as an observable — analysis_overview §3, roadmap §Q2.9 (in-scope is an **open** question; machinery "not assumed") | [method-we-use] (conditional / open scope) |
| **ρ = (p_ID − p_MS)/p_ID momentum-imbalance template fit** to separate heavy-flavour muons from π/K decay-in-flight + punch-through ("fake") muons | Our fake-muon **purity** procedure — the Δp/p significance template fit (roadmap task_07, IntNote §16 "Background / fake-muon purity"); analysis_overview §6 states this method is "adapted from the Run 2 note." This single-muon paper is the origin of that method | [method-we-use] |
| Per-muon efficiency weighting w = 1/(ε_reco·ε_trig), tag-and-probe J/ψ→μμ, binned in p_T and **charge×η (q·η)** | Background/comparison for our trigger-efficiency write-up (we also bin single-muon ε^nc in q·η) — IntNote §8/§10; **but** our reco correction is a PAIR efficiency (see warning) | [background-for-writing] |
| d_0 (transverse impact parameter) template fit separating **charm vs bottom** muons | Intro / thesis background on c-vs-b separation; NOT our signal definition (we select a single-b dimuon) | [background-for-writing] |
| Physics motivation: b-quark flow / dead-cone mass hierarchy (charm v_n > bottom v_n) as a QGP probe | Introduction / thesis ch.1 — motivates b-hadron flow as a potential observable | [background-for-writing] |
| Measured inclusive/charm/bottom muon v_2, v_3 vs p_T and centrality | A prior measurement to cite/compare if we measure dimuon flow | [background-for-writing] |

## Scope & condition-difference warnings

Pb+Pb **5.02 TeV**, ATLAS Run 2, 2015 (0.49 nb⁻¹) + 2018 (1.44 nb⁻¹) = up to
1.9 nb⁻¹ [PLB p.1,3]. Signal = **inclusive single heavy-flavour muons**
(charm OR bottom), 4 < p_T < 30 GeV, |η| < 2.0.

- **Different signal object.** This measures *single* HF muons (one muon per
  heavy quark, charm and bottom from *two different* sources). OUR signal is a
  **dimuon from a single b-hadron** (b→μ + b→c→μ, low mass, small opening
  angle). Do not conflate the selections or the c/b separation.
- **Single-muon ≠ pair efficiency.** Reco efficiency here is factorized
  single-muon ε_reco = ε_ID(data)·S_MS·ε_MS^MC [note eq.15]. OUR reco
  correction is a **3-D pair efficiency** ε_reco(pair p_T, pair η, ΔR), NOT
  ε₁·ε₂ (analysis_overview §4b; project memory `pair_reco_effcy_definition`).
  Do not import the single-muon factorization for our reco step.
- **v_n is an OPEN observable for us.** Our primary results are cross-sections
  and R_AA; flow (v_2) is a *potential* addition only (overview §3, roadmap
  §Q2.9). Treat the EP machinery here as a reference design, not a committed step.
- **Run 2 5.02 TeV vs our Run 3 5.36 TeV / HIJING overlay.** Detector
  conditions, occupancy, and trigger menu differ. ACKNOWLEDGE; do NOT assume
  the size or direction of any performance difference (GUIDE §5).

## Content summary

### HF-muon selection & efficiency  [PLB p.3; note §3.4, §6.2]
- Muons reconstructed in both ID and MS, **'medium' WP**, 4 < p_T < 30 GeV,
  |η| < 2.0, **no TRT-hit requirement** (high HI occupancy). Each muon matched
  to a trigger muon; weight = inverse of (reco × trigger) efficiency.
- Triggers: single-muon HLT thresholds p_T > 3, 4, 6, 8 GeV; p_T>8 unprescaled
  (full lumi), lower thresholds prescaled. Sampled lumi 0.3 / 0.6 / 1.9 nb⁻¹
  for muons in 4–6 / 6–8 / >8 GeV.
- Reco eff [note eq.15]: ε_reco = ε_ID(data, J/ψ T&P) · S_MS^reco ·
  ε_MS^MC(Medium). ID eff ≈ 99% all centralities (from Pb+Pb data); MS eff from
  J/ψ data-overlay MC with data/MC scale factors. Maps in p_T vs **q·η**
  (charge×η used because low-p_T reco is charge-dependent).
- Trigger eff [note eq.16]: ε_trig = ε_trig^MC(signal-only Υ(nS) Pythia8) ·
  S_trig^pp(2017 pp data/MC SF) · f_cent (centrality-dependent factor from
  Pb+Pb/pp comparison). Data/MC corrections ~1–10%.

### Centrality  [PLB p.3]
FCal ΣE_T percentiles: 0–10, 10–20, 20–30, 30–40, 40–60%. MC Glauber ⟨N_part⟩ =
358.8±2.3, 264.1±2.9, 189.2±2.8, 131.4±2.6, 70.5±2.2 respectively.

### Event-plane method (flow extraction)  [PLB p.4; note §5, §6.4]
- FCal q-vector per event [note eq.3]: q_n = Σ_i E_T(φ_i)(cos nφ_i, sin nφ_i) /
  ΣE_T; EP angle Ψ_n = (1/n) arctan(q_{n,y}/q_{n,x}) [eq.4]. q recentered
  (subtract ⟨q⟩) and flattened via inverse-sqrt covariance matrix [eqs.5–7],
  per FCal side, calibrated on minimum-bias events.
- **Resolution** Res{nΨ_n} by two-subevent method (FCal A vs C) [note eq.8]:
  Res{nΨ_n^{A|C}} = √⟨cos n(Ψ_n^A − Ψ_n^C)⟩; full-FCal resolution from
  χ_n^Full = √2·χ_n via the Bessel-function relation [eq.9]. Res rises with
  centrality, consistent between 2015 and 2018.
- **Raw v_n**: yield N_X^μ fit vs n(φ−Ψ_n) to (1/N)dN/d(n(φ−Ψ_n)) =
  1 + 2 v_n^raw cos(n(φ−Ψ_n)) [PLB eq.; note eq.17], X ∈ {inclusive, charm,
  bottom}, n = 2, 3. **Corrected** v_n = v_n^raw / Res{nΨ_n} [eq.18].

### Signal/background separation — two-step ρ then d_0 template fit  [PLB p.5; note §6.1–6.3]

> **Shared mechanism (ρ=Δp/p; d_0) lives in the concept hub
> [[../concepts/muon_source_template_fits]].** Below = *this paper's specifics*
> (2D PDF factorization, template provenance, shift/smear, η intervals, systematics).

Two discriminating variables:
- **ρ = (p_ID − p_MS)/p_ID** (momentum imbalance; p_MS energy-loss-corrected).
  Real muons peak at ρ≈0; π/K background is broader, shifted positive.
- **d_0** transverse impact parameter wrt primary vertex; charm vs bottom differ
  (decay length).

Backgrounds: (1) **π/K** decay-in-flight + punch-through; (2) **light/onia**
(direct quarkonia, low-mass resonances, τ, W/Z) — high-quality muons subtracted
from signal; (3) **Calo** (hadronic showers, mismatched/punch-through) — small.

2D PDF = Σ_{i=1..5} κ_i P_i(ρ)⊗D_i(d_0) [note eq.10]; ρ and d_0 found
uncorrelated for signal → factorized into a **two-step 1-D fit** [eqs.11–13]:
1. **ρ fit** → inclusive HF-muon yield (N_c+N_b) and π/K yield. Light/onia
   fraction **fixed** to Pythia8 (~4% avg); Calo ignored.
2. **d_0 fit** (with π/K fixed) → split inclusive into charm N_c and bottom N_b.

Template provenance [PLB p.5; note §6.3]:
- Signal (charm/bottom) ρ & d_0 from Pythia8 multijet hard-scattering (A14 tune,
  NNPDF2.3LO) filtered on a generator-level muon; charm and bottom ρ shapes
  identical.
- π/K templates from non-diffractive QCD Pythia8; light/onia from direct J/ψ.
- **ρ templates: no Pb+Pb overlay** (overlay MC ID resolution worse than data);
  single-muon ID/MS momentum response **shifted & smeared** to match data, tuned
  on J/ψ→μμ invariant-mass resolution vs centrality.
- **d_0 templates: with MB Pb+Pb overlay** (to reproduce vertex resolution);
  shifted/smeared on high-quality prompt-track d_0; charm/bottom signal
  reweighted to measured (ALICE/CMS) c/b-hadron p_T spectra (R_AA-based).
- Fits done independently in p_T, centrality, n|φ−Ψ_n|, and **two η intervals**
  (|η|<1, 1<|η|<2), then combined. Likelihood within |d_0|<0.5 mm (~97% of
  muons), normalization within |d_0|<1.0 mm [note §6.3.3].

### Systematics  [PLB p.6–9, Table 1]
Eight categories: (1) muon efficiency, (2) ρ fit, (3) d_0 fit, (4) light/onia
bkg, (5) other bkg, (6) ρ–d_0 correlation, (7) EP resolution, (8) jet bias.
EP-resolution and jet-bias found negligible, not in final totals; summed in
quadrature. Typical **absolute** sizes for v_2 (v_3) [Table 1]:

| Category | Inclusive HF v_2(v_3) | Charm v_2(v_3) | Bottom v_2(v_3) |
|---|---|---|---|
| Muon efficiency | <0.0002 (0.0006) | <0.006 (0.001) | <0.001 (0.001) |
| ρ fit | <0.004 (0.006) | <0.009 (0.01) | <0.005 (0.003) |
| d_0 fit | N/A | <0.02 (0.03) | <0.01 (0.01) |
| Light/onia bkg | <0.004 (0.002) | <0.02 (0.01) | <0.008 (0.004) |
| Other bkg | <1e-6 (1e-6) | <0.002 (0.004) | <0.001 (0.0004) |
| ρ–d_0 correlation | N/A | <0.01 (0.004) | <0.007 (0.005) |

- ρ-fit syst: vary shift/smear (none; signal-only; extra bkg smearing in ρ>0.2);
  1–10%, largest at low p_T / high η, no centrality dependence.
- d_0-fit syst: 2018 vs 2015 shift/smear params + R_AA reweighting variation;
  1–20%, largest at high p_T / peripheral. d_0-fit systs **anti-correlated**
  between charm and bottom.
- Light/onia: vary fixed fraction (0 and 2× Pythia) and its assumed v_n (0 or 2×
  inclusive); anti-correlated charm↔bottom.
- ρ–d_0 correlation: redo d_0 fit with ρ<0.1 cut.
- Jet bias (recoil jet aligning Ψ_n): 0.3–0.4%, negligible.

### Results  [PLB p.10–14]
- Inclusive HF-muon v_2, v_3 measured 4 < p_T < 30 GeV in 5 centrality bins;
  both **decrease with p_T**. v_2 smaller in central (0–10%) as expected from
  smaller initial ellipticity; v_3 nearly centrality-independent (fluctuation
  origin). Consistent with prior ATLAS 2.76 TeV results.
- Charm and bottom v_2, v_3 separated **only up to p_T = 20 GeV** (above that the
  c/b separation is unstable). **Charm v_2 (v_3) > bottom v_2 (v_3)** at all p_T
  and centrality — qualitatively matching mass/dead-cone expectations (heavier b
  less modified). Charm/bottom uncertainties partially **anti-correlated**.
- Compared to DREENA-B and DAB-MOD heavy-quark energy-loss/transport models.

## References worth future reading   (§6; ≤3)
1. **ATLAS, PRC 98 (2018) 044905, arXiv:1805.05220** [ref 7] — *Suppression and
   azimuthal anisotropy of HF-decay muons in Pb+Pb 2.76 TeV*. PRIMARY. New info:
   the **R_AA** measurement for HF muons (directly our primary observable's
   single-muon analogue) plus the original event-plane method this paper extends.
   Serves R_AA (overview §3) + flow. NOTE: a companion R_AA **internal note** PDF
   is already on disk (`Run2 HF muons nuclear modification factor internal
   notes.pdf`, not yet summarized) — coordinate to avoid duplication.
2. **ATLAS, PRL 124 (2020) 082301, arXiv:1909.01650** [ref 33] — *Azimuthal
   anisotropy of charm/bottom muons in pp 13 TeV*. SUPPORTIVE. New info: the **pp
   baseline** for HF-muon flow and the same ρ/d_0 template method applied in pp;
   relevant if we compare dimuon flow to a pp reference. Serves Intro/method.
- **Not added:** ATLAS muon reco performance, EPJC 76 (2016) 292,
  arXiv:1603.05598 [ref 40] (the tag-and-probe ε method used here) — already
  covered by the newer KB doc `physics/detector/ATLAS_Run2_muon_reconstruction.md`
  (EPJC 81 (2021) 578).

## Related KB docs   (knowledge graph)
- [[../concepts/muon_source_template_fits]] — the ρ=Δp/p (real vs π/K) and d0
  (charm vs bottom) template fits this paper uses, as a cross-paper concept hub.
- [[run2_hf_muon_raa]] — companion R_AA measurement on the same muon sample.
- [[yun_tian_thesis]] — PhD thesis with the in-depth (2.76 TeV) method companion (EP + scalar-product v_n) to this HF-muon v_n analysis.
- [[physics/heavy_ion/open_hf_production]] — field background (energy loss, flow, dead-cone).
- [[run2_dimuon_note.md]] — Run 2 *dimuon* note (ATL-COM-PHYS-2021-1094); uses
  the **same ρ=Δp/p momentum-imbalance template fit** for fake-muon purity
  (pair-significance version) and FCal centrality; that note is the direct
  predecessor of OUR analysis, this paper is the single-muon origin of the
  Δp/p method.
- [[physics/detector/ATLAS_Run2_muon_reconstruction.md]] — ATLAS muon
  reco/ID/trigger tag-and-probe efficiency methodology (newer than ref 40 here).
- [[analysis/overview.md]] — physics goal, where v_2 sits as an open observable.
</content>
</invoke>
