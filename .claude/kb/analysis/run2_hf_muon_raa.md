# Run 2 Heavy-Flavour Muon R_AA (charm/bottom separated)

**Source:** ATLAS Collaboration, *Measurement of the nuclear modification
factor for muons from charm and bottom hadrons in Pb+Pb collisions at
5.02 TeV with the ATLAS detector* — Phys. Lett. B 829 (2022) 137077
(CERN-EP-2021-153).
**arXiv / DOI:** arXiv:2109.00411 ; 10.1016/j.physletb.2022.137077
**PDF:** `./2109.00411.pdf` (PRIMARY, journal Letter)
**Companion internal note (read for implementation detail):**
`./Run2 HF muons nuclear modification factor internal notes.pdf`
= ANA-HION-2019-58-INT1 (Q. Hu, D. Perepelitsa, J. Nagle), 4 Apr 2021,
v1.0. The note holds the full method (selections, momentum-imbalance and
d0 template fits, efficiencies, systematics); the Letter holds the
results and the condensed method.
**Classification:** PRIMARY
**Added:** 2026-06-14

---

## Relevance to this analysis   (this is the paper OUR analysis is built upon)

OUR analysis (Run 3 single-b → **dimuon**) directly descends from this
single inclusive-HF-muon measurement. The two share the physics question
(b-quark in-medium energy loss via R_AA), much of the detector method
(muon selection, efficiency corrections, FCal/Glauber centrality, the
momentum-imbalance fake-muon separation), and the physics motivation
(dead-cone, mass-ordered energy loss). OUR value-add: (1) **dimuon pair
kinematics** (pair pT, pair η) give better kinematic control of the parent
b than one HF muon; (2) the **single-b dimuon** topology (low mass, small
opening angle) isolates b-decay products from c-hadron contamination →
flavour-pure signal *by selection*, instead of via this paper's per-muon
d0 template fit.

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| **R_AA definition + per-event-yield normalization** (Eq. 1, 3): R_AA = (N_AA/N_evt) / (⟨T_AA⟩·σ_pp) | R_AA observable — `analysis_overview` §3, roadmap step 7 (`task_06_raa`), IntNote §15 | [method-we-use] |
| **FCal ΣE_T + MC-Glauber centrality**, ⟨T_AA⟩/⟨N_part⟩ per bin | Centrality determination — overview §4e, IntNote §6 (method; our numbers are 5.36 TeV, differ) | [method-we-use] (method) / [background-for-writing] (their numbers) |
| **Momentum-imbalance ρ=(p_ID−p_MS)/p_ID template fit** to separate real (HF) muons from hadronic/fake background | Δp/p significance fake-muon purity framework — overview §6, roadmap step 8 (`task_07`), IntNote §16 | [method-we-use] |
| **Per-(single)muon efficiency weight** w=1/(ε_reco·ε_trig), tag-and-probe; reco-eff = ε(Medium\|ID)·ε(ID\|ME) factorization | Efficiency-correction framing — IntNote §8–12 (our reco-eff is a **pair** eff, see warnings) | [background-for-writing] (method reference) |
| **HF single-muon selection** (Medium WP, 4<pT<30 GeV, \|η\|<2, trigger-matched) | Single-muon selection / cutflow — IntNote §7 | [background-for-writing] |
| **d0 (impact-parameter) template fit for charm↔bottom separation** (B vs D lifetimes) | Flavour-separation methodology to cite/contrast (we separate via dimuon kinematics) — Intro / thesis | [background-for-writing] |
| **Prompt-muon background** (quarkonia, W/Z) sizes & estimation | Background discussion for the signal region — Intro / IntNote §16 | [background-for-writing] |
| **Physics motivation**: dead-cone, mass-ordered (c more suppressed than b) energy loss; R_AA+v2 joint constraint | Introduction / thesis ch. 1 | [background-for-writing] |
| **pp cross-section vs FONLL**; charm/bottom muon spectra | Precedent for pp24 crossx-vs-theory (we use POWHEG+Pythia) — IntNote §14 | [background-for-writing] |
| **Systematics source list** (muon eff / background removal / c-b separation / global norm) | Systematics framework — roadmap `task_08`, IntNote §17 | [background-for-writing] |

---

## Scope & condition-difference warnings

System/energy/dataset: **Run 2, √s_NN = 5.02 TeV.** pp = 2017, 1.17 pb⁻¹.
Pb+Pb = 2015 (208 µb⁻¹) + 2018 (38 µb⁻¹) = **246 µb⁻¹ combined**.

Differences from OUR analysis — **acknowledge; do NOT assume size or
direction (GUIDE §5):**
- **Run 2 5.02 TeV** here vs **Run 3 5.36 TeV** (pp24 + Pb+Pb 23/24/25) ours.
- **Single inclusive HF muon** here vs **single-b dimuon pair** ours.
  Their flavour separation is a per-muon d0 template fit; ours is a
  kinematic selection on the pair (low mass, small ΔR).
- **Trigger:** single **HLT_mu4** (L1+HLT pT>4 GeV) for both pp and Pb+Pb
  here; ours is PbPb single-mu4 and **pp24 2mu4** (different efficiency
  machinery).
- **Reconstruction efficiency:** here a **single-muon** efficiency,
  factorized ε(Medium|ID)·ε(ID|ME) from J/ψ/Z tag-and-probe. OURS is a
  **pair** efficiency in (pair pT, pair η, ΔR) from Pythia fullsim /
  HIJING overlay — *not* a product of single-muon efficiencies (see
  `project_pair_reco_effcy_definition`). These are fundamentally
  different objects.
- **Muon WP:** Medium (nominal) here; our nominal WP is unsettled
  (roadmap Q2.8).
- **Centrality:** Glauber with σ_NN = 70 mb at 5.02 TeV here; our PbPb
  uses σ_PbPb = 7.8 b at 5.36 TeV (roadmap Q2.2). T_AA values below are
  5.02 TeV and must NOT be reused for Run 3.

---

## Content summary

### Physics motivation  [background-for-writing]
Heavy quarks (c, b) are made in initial hard scatters, conserve flavour
(cannot be destroyed in the QGP, only have momenta modified), and probe
the QGP via radiative + collisional energy loss whose balance depends on
quark mass. Radiative loss is suppressed by the **dead-cone effect**
(gluon radiation suppressed at angles < m/E), more so for b than c at
equal momentum. Measuring c and b R_AA separately constrains the
mass-dependence of energy loss (§1). Joint R_AA + v2 (the latter from the
companion measurement, ref [25]) discriminates energy-loss mechanisms and
QGP-geometry effects.

### R_AA and yields  [method-we-use]
R_AA (Eq. 1):  **R_AA = (N_AA / N_evt) / (⟨T_AA⟩ × σ_pp)** , where N_AA =
observed muons in Pb+Pb, N_evt = number of min-bias Pb+Pb events,
⟨T_AA⟩ = mean nuclear thickness, σ_pp = pp production cross-section at the
same energy. R_AA = 1 ⇒ scaled pp; R_AA < 1 ⇒ suppression.
pp differential cross-section (Eq. 2): d²σ/dpT dη = N_corr /
(ΔpT·Δη·L_HLT_mu4). Pb+Pb per-event yield (Eq. 3): N_AA =
N_corr / (ΔpT·Δη·(N_evt^MB/L_MinBias)·L_HLT_mu4), with the trigger and
min-bias sampled luminosities cancelling. N_corr = background-subtracted
muon counts weighted per muon by **w = 1/(ε_reco·ε_trig)**.
(Internal note §4.1; N_evt^MB/L_MinBias = 7.383×10⁸ nb⁻¹ per 10%.)

### Event / muon selection  [background-for-writing]
Good event: detector-error flags clear, ≥1 primary vertex; Pb+Pb only:
\|Δt_MBTS\|<5 ns and FCal–ZDC / calo-Σtransverse-energy–nTrk pileup
rejection (rejects ~0.2%); HLT_mu4 fired. Muon: **Medium** WP
(TRT-hit requirement in pp only), **4 < pT < 30 GeV, \|η\| < 2**,
\|z0 sinθ\| < 1 mm, matched to online HLT_mu4 muon with ΔR < 0.01
(note §3.3).

### Centrality  [method-we-use (method)]
Per-event ΣE_T^FCal percentile, MC Glauber (σ_NN = 70 mb). Bins and
5.02 TeV Glauber values (note Table 5; T_AA in mb⁻¹):

| Centrality | ⟨N_part⟩ | ⟨T_AA⟩ [mb⁻¹] |
|---|---|---|
| 0–10%  | 358.8 | 23.35 |
| 10–20% | 264.1 | 14.33 |
| 20–30% | 189.2 | 8.64 |
| 30–40% | 131.4 | 4.95 |
| 40–60% | 70.5  | 1.96 |

### Background composition & two-step template fit  [method-we-use (ρ fit); background (d0)]

> **Shared mechanism (ρ=Δp/p real-vs-fake; d0 charm-vs-bottom) lives in the
> concept hub [[../concepts/muon_source_template_fits]].** Below = *this paper's
> specifics* (component model, template provenance, FONLL reweighting, fit ranges).

Selected muons (4<pT<30, \|η\|<2) = signal HF muons (c+b semileptonic) +
three backgrounds: **prompt-muon** (real muons from quarkonia, W/Z, τ —
fixed from simulation), **hadronic** (π/K decay-in-flight or punch-through
in the ID/calo), **fake** (random ID+MS segment combinations; fixed
relative to hadronic). Two discriminating variables, found uncorrelated
for signal, so the 2-D PDF factorizes (note Eq. 7–12):

1. **Momentum-imbalance fit.** ρ = (p_ID − p_MS)/p_ID (p_MS energy-loss
   corrected). Real (HF + prompt) muons: symmetric peak at ρ≈0; hadronic
   background: broad, shifted to positive ρ (energy loss not properly
   modelled). 4-component template fit (signal P_sig(ρ), hadronic
   P_had(ρ), fake P_fake(ρ), prompt fixed) — effectively 2 free
   components, since prompt and fake/hadronic ratio are externally fixed.
   Output: HF-muon yield separated from hadronic+fake. Templates from
   muon-filtered/non-diffractive **Pythia8** QCD at 5.02 TeV (A14,
   NNPDF23lo); muon momentum calibrated to data via shift Δp / smear δp
   factors from J/ψ→µµ invariant-mass response (note §4.4.1).
2. **d0 (impact-parameter) fit.** Within the fixed HF yield, charm vs
   bottom separated using the muon track transverse impact parameter d0
   (relative to the **beam spot**) — different B/D lifetimes (B≈1.5×10⁻¹²,
   D⁺≈1.0×10⁻¹², D⁰≈0.4×10⁻¹² s) give different d0 spreads. Two-component
   fit over **\|d0\| < 0.5 mm**; the single free parameter is the
   **bottom fraction f_b = κ_bottom/(κ_charm+κ_bottom)**. Bottom muons
   include the b→c→µ cascade. Signal d0 templates: muon-filtered Pythia8,
   reweighted to **FONLL** c/b-meson pT spectra, with baryon-to-meson and
   D⁺/D⁰ corrections; in Pb+Pb additionally reweighted to ALICE/CMS
   measured modified hadron spectra and built by overlaying Pythia8 onto
   2015 min-bias Pb+Pb events.

Prompt-muon background size: dominated by W decays at high pT — ~3% in pp
and ~10% (5%) in 0–10% (40–60%) Pb+Pb at 20<pT<30 GeV (Letter §4).

### Efficiency corrections  [background-for-writing (method reference)]
- **Reconstruction** (note §4.6): single-muon ε_reco = ε(Medium|ID) [MS
  reco+ID, "MS efficiency"] × ε(ID|ME) [ID efficiency], measured by
  J/ψ→µµ / Z tag-and-probe. pp efficiency from large **13 TeV** Muon-CP
  T&P maps, applied to 5.02 TeV (difference negligible, ~0.3%); MS eff
  plateaus ~97% at pT>7 GeV. Pb+Pb ID eff ~98% from Pb+Pb J/ψ T&P; Pb+Pb
  reco eff factorized as ε_pp × ε_PbPb(ID|ME)/ε_pp(ID|ME) (pp ID eff
  ≈100% so dropped). **No centrality dependence** observed for MS or ID
  reco efficiency.
- **Trigger** (note §4.5): central value from J/ψ→µµ (and Υ(nS)) Pythia8
  simulation tag-and-probe vs (pT, q·η); data/MC mismodelling corrected
  by a scale factor ≤10%. Pb+Pb adds a data-driven A_trig term (Pb+Pb/pp
  efficiency ratio) for online-momentum-scale and 2015↔2018 differences,
  as a function of pT, η, centrality, year.

### Systematic uncertainties  (Letter Table 1; ranges over pT & centrality)  [background-for-writing]

| Source | σ_pp c→µ [%] | σ_pp b→µ [%] | R_AA c→µ [%] | R_AA b→µ [%] |
|---|---|---|---|---|
| Muon efficiency | 0.5–1.0 | 0.4–0.6 | 0.3–16 | 0.2–16 |
| Background removal | 4.3–12 | 0.8–3.8 | 1.9–27 | 0.5–4.7 |
| Charm–bottom separation | 4.5–9.8 | 3.2–8.0 | 5.1–23 | 4.1–13 |
| Global normalization | 1.6 | 1.6 | 1.8–4.9 | 1.8–4.9 |
| **Total** | **5.8–13** | **3.1–7.4** | **6.5–35** | **5.4–18** |

Leading sources: background removal and charm–bottom separation
everywhere; muon efficiency also leads at low pT in central Pb+Pb. pp
luminosity uncertainty 1.6% (LUCID-2); ⟨T_AA⟩ from Glauber-parameter
variation. Some sources cancel partially in R_AA (pp↔Pb+Pb correlated).

### Results  [background-for-writing]
- **pp cross-sections vs FONLL:** b→µ (incl. b→c→µ) is 1.3–1.4× the FONLL
  central value but within its band at low pT, agreeing for pT>10 GeV;
  c→µ lies above FONLL at all pT, consistent within FONLL uncertainty
  below 10 GeV.
- **R_AA:** significant suppression of both c→µ and b→µ in all
  centralities; suppression increases monotonically peripheral→central.
  c→µ more suppressed than b→µ at low pT (<10 GeV) — mass-ordering
  consistent with dead-cone expectation. Charm R_AA ≈ flat in pT; bottom
  R_AA dips to ~10 GeV then flattens. Data compared with DREENA-B and
  DAB-MOD energy-loss models; joint with v2 from ref [25].

---

## References worth future reading   (§6; ≤3)

1. **ATLAS, *Measurement of azimuthal anisotropy of muons from charm and
   bottom hadrons in pp collisions at √s=13 TeV*, PRL 124 (2020) 082301,
   arXiv:1909.01650** — *PRIMARY candidate.* Origin of the **d0
   charm/bottom template-separation** method in pp; the detailed
   reference for the impact-parameter approach this paper reuses. Serves:
   flavour-separation methodology + Intro.
2. **ATLAS, *Measurement of azimuthal anisotropy of muons from charm and
   bottom hadrons in Pb+Pb collisions at √s_NN=5.02 TeV*, PLB 807 (2020)
   135595, arXiv:2003.03565** — *PRIMARY candidate.* Companion **v2**
   measurement on the same muon sample; directly relevant to our
   *potential* dimuon v2 observable (overview §3.3) and the source of the
   muon momentum shift/smear calibration. Serves: v2 observable + flow
   background.
3. **ATLAS, *Suppression and azimuthal anisotropy of HF muons in Pb+Pb at
   √s_NN=2.76 TeV*, PRC 98 (2018) 044905, arXiv:1805.05220** —
   *SUPPORTIVE.* Predecessor inclusive-HF-muon R_AA (no c/b split);
   background-for-writing only (history/Intro).

---

## Related KB docs

- [[../concepts/muon_source_template_fits]] — the Δp/p (real vs π/K) and d0
  (charm vs bottom) template fits this paper uses, as a cross-paper concept hub.
- [[run2_hf_muon_vn]] — companion v_n measurement on the same muon sample.
- [[physics/heavy_ion/open_hf_production]] — field background (R_AA, dead-cone, energy loss).
- [[atlas_centrality_2023]] — source of our ⟨T_AA⟩ values (5.36 TeV); [[glauber_modeling]] — the underlying centrality methodology this paper's R_AA normalization uses.
- [[run2_dimuon_note]] — the Run 2 **dimuon** Δφ analysis from the same
  group; shares the ρ=Δp/p momentum-imbalance fake-muon template fit and
  the FCal/Glauber centrality. (Different observable: away-side widths.)
- [[overview]] — our analysis goal, signal, R_AA, sample roles.
- [[physics/detector/ATLAS_Run2_muon_reconstruction]] — the tag-and-probe
  muon reco/trigger efficiency methodology this paper applies.
</content>
</invoke>
