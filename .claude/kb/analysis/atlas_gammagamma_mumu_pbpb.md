# ATLAS γγ→μμ in non-UPC Pb+Pb (5.02 TeV) — per-pair trigger/reco efficiency method

**Source:** ATLAS Collaboration, *Measurement of muon pairs produced via γγ
scattering in non-ultraperipheral Pb+Pb collisions at √s_NN = 5.02 TeV with the
ATLAS detector* — Phys. Rev. C **107** (2023) 054907; CERN-EP-2022-047.
**arXiv / DOI:** 2206.12594 [nucl-ex] / 10.1103/PhysRevC.107.054907
**PDF:** `./2206.12594.pdf` (PRIMARY)
**Classification:** PRIMARY
**Added:** 2026-06-15

> This is **Ref. [60]** of the Run-2 dimuon back-to-back Letter
> ([[run2_dimuon_backtoback_paper]]) — it is the **method provenance** for the
> "per-pair inverse trigger×reco efficiency weight" that the Letter only *cites*.
> The Letter's machinery (single-b away-side correlation) is summarized there; the
> per-muon efficiency definitions, the trigger combination, and the γγ→μμ
> identification cuts live **here**. Use both; do not duplicate.
>
> **CRITICAL contrast — read before reusing anything:** this paper's dimuon pairs
> are **back-to-back and well-separated in the detector**, which is the explicit
> physical justification for factorizing the pair reco efficiency as a **product
> of two single-muon efficiencies** ε_rec1·ε_rec2 (and for combining single-muon
> trigger efficiencies). **OUR signal is the opposite** — a single-b pair of
> **nearby, close-by muons** that are reconstruction-correlated — so our reco
> efficiency is a single **PAIR** object ε_reco(pair pT, pair η, ΔR), **not**
> ε₁·ε₂ (analysis_overview §4b). The trigger *combination* method transfers (we
> add a ΔR correlation factor); the **reco-eff factorization does NOT**.

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Per-pair trigger efficiency built by **combining independent single-muon trigger efficiencies** ε(pT, q·η), measured by matching offline muons to a trigger object; done separately per trigger chain (L1Single / L1Pair) — §4.1 | Trigger-eff method provenance for roadmap **step 8** (PbPb single-mu4) / **step 10** (pp24 2mu4); IntNote §8. We transfer the single-μ→pair combination but **add a ΔR correction** (overview §4b) | [method-we-use] |
| Pair **reconstruction efficiency = product of single-muon reco efficiencies** ε_rec1·ε_rec2, fn of (pT, q·η, centrality), from STARlight+HIJING overlay MC, justified by back-to-back / well-separated topology — §4.1 | The actual definition the back-to-back Letter cites. **CONTRAST for roadmap step 12:** our pairs are nearby → we use a 3D PAIR ε_reco(pair pT, pair η, ΔR), NOT ε₁·ε₂ (overview §4b) | [method-we-use] |
| Per-pair efficiency weight w = 1/(ε_trig·ε_rec1·ε_rec2·ε_vtx) applied event-by-event; plus a **vertex-pointing efficiency ε_vtx** for the d0pair / z0sinθpair cuts — §4.1 | Per-pair weighting structure (overview §4b); our weight is w⁻¹ = ε_trig^pair·ε_reco^pair (we have no separate vertex-pointing cut) | [method-we-use] |
| **γγ→μμ identification** via pair **acoplanarity α** and **asymmetry A** (and the k⊥ scale) — §1, §3.2 | γγ→μμ is a **background we REMOVE** from OS pairs ([[run2_dimuon_note]]); this paper *selects* the γγ region that we *veto*. Pair-selection definitions (IntNote §7) | [method-we-use] |
| HF-decay = the **dominant background** here, suppressed by muon-pair vertex pointing (d0pair, z0sinθpair) + a **d0pair template fit** (MC signal template, data sideband background template) — §3.2, §4.2 | Their HF background is **our signal**. Complementary to our Δp/p purity ([[muon_source_template_fits]], roadmap step 16); shows a d0-based separation alternative | [background-for-writing] |
| DY (OS prompt-dimuon) background via Powheg+Pythia8 with 5 nuclear PDF sets (nNNPDF2.0 nominal), normalized via effective NN cross-section × ⟨T_AA⟩ — §3.1, §4.4 | OS DY treatment if/when needed; IntNote backgrounds | [background-for-writing] |
| FCal ΣE_T Glauber centrality, **⟨T_AA⟩ table at 5.02 TeV** — §3.3, Table 1 | Centrality method (overview §4e) and a 5.02-TeV T_AA reference (our values are 5.36 TeV — differ) | [background-for-writing] |
| γγ→μμ as a coherent-photon **QED process / EM probe of the QGP** (acoplanarity broadening, magnetic-field tests) | Intro / thesis context — a non-UPC QED dilepton process; explicitly **not** our HF signal | [background-for-writing] |

## Scope & condition-difference warnings

Run 2, √s_NN = **5.02 TeV** Pb+Pb (2015 + 2018), **1.93–1.94 nb⁻¹**. Performance
MC = **STARlight γγ→μμ overlaid on HIJING v1.383 min-bias**, Geant4 (4M events).
Muons use the **"medium" working point** (Ref. [41]) as nominal — note this
**differs from the back-to-back Letter, which used "tight"**; tight is here only a
systematic (yields drop ~20%). Single muon: pT > 3.7 GeV, |η| < 2.4, m_μμ < 45 GeV.

Our analysis is **Run 3, 5.36 TeV** (pp24 + PbPb 23/24/25). ACKNOWLEDGE the
energy / conditions / detector / pileup difference; **DO NOT** estimate its size
or direction (GUIDE §5). The T_AA values below are 5.02 TeV — our analysis uses
5.36 TeV Glauber numbers (roadmap §Q2.3); do not transfer them.

**Topology difference (the load-bearing one):** the factorized reco efficiency
(ε₁·ε₂) and the single-μ trigger combination are valid **because these pairs are
back-to-back / well-separated**. Our nearby single-b muons are
reconstruction-correlated → pair ε_reco(pair pT, pair η, ΔR). Do not infer that
ε₁·ε₂ is acceptable for our signal.

## Content summary

**Physics (Intro, §1) [background-for-writing].** Lorentz-contracted EM fields of
the Pb ions act as a source of quasi-real photons (Weizsäcker–Williams /
equivalent-photon approx); γγ→μ⁺μ⁻ occurs even in **non-UPC** (hadronic-overlap)
collisions. The dimuon acoplanarity/k⊥ broadening with centrality is studied as a
possible EM probe of the QGP (rescattering off plasma constituents, long-range
magnetic fields). This is a **QED process, not HF** — the inverse of our analysis,
where HF dimuons are the signal and γγ is a contaminant.

**Kinematic variables (§3.2).** With φ₁,₂ and pT1,2 the muon azimuths and momenta:
- **Acoplanarity** α ≡ 1 − |φ₁ − φ₂|/π   (→ 0 for back-to-back pairs)
- **Asymmetry** A ≡ |pT1 − pT2| / (pT1 + pT2)
- **k⊥ (transverse momentum kick)** ≡ ½(pT1+pT2)(π − |φ₁−φ₂|) = α·π·p̄T,
  with **p̄T ≡ (pT1+pT2)/2** (Eq. 1). γγ pairs have very small α, A, k⊥.

**γγ→μμ fiducial selections (§3.2).** Two regions (each also includes the η/pT/mass
preselections):
- **Fid-α**: A < 0.06 ∧ α < 0.012   (used for acoplanarity distributions)
- **Fid-k⊥**: A < 0.06 ∧ k⊥ < 150 MeV   (used for k⊥ distributions and the
  cross-sections)
Separate regions are used because α and k⊥ are directly related (Eq. 1), so a
single cut would distort the p̄T-dependence of the HF background. **Note for our
analysis:** the back-to-back Letter / [[run2_dimuon_note]] *veto* this γγ region
(keep pairs with α≳0.01 OR A≳0.08); this paper *selects* it.

**Vertex pointing (§3.2).** To kill HF (displaced) muons:
d0pair ≡ √(d0₁² + d0₂²) < 0.1 mm and (z0 sinθ)pair ≡ √((z0sinθ)₁² + (z0sinθ)₂²)
< 0.2 mm. Reduces HF pairs by ≈2× at a γγ-signal inefficiency < 2%.

**Per-pair efficiency correction (§4.1) — the key method.** Each data pair is
weighted by
```
w = 1 / ( ε_trig · ε_rec1 · ε_rec2 · ε_vtx )
```
- **ε_trig (pair trigger eff):** because the pairs are **back-to-back and
  well-separated**, the dimuon-trigger efficiency is built from **independent
  single-muon trigger efficiencies**, measured by testing whether an
  offline-reco muon (passing preselection) is **matched to a trigger-found muon**,
  as a function of muon **pT and q·η** (charge × pseudorapidity). The single-muon
  efficiencies are combined **separately for the two dimuon triggers** — L1Single
  (1μ pT>4 @L1 + 2μ pT>4 @HLT) and L1Pair (2μ pT>4 @L1 and @HLT) — into per-pair
  trigger efficiencies. Trigger eff varies only **a few %** across the full
  centrality range (occupancy varies by orders of magnitude).
- **ε_rec1·ε_rec2 (pair reco eff):** taken as the **product of single-muon reco
  efficiencies**, from the STARlight+HIJING overlay MC, fn of muon **pT and q·η**,
  corrected for small data–MC differences (Ref. [41]). **Negligible** centrality
  dependence at mid-rapidity (|q·η| < 1); **≈10% decrease** peripheral→central at
  forward |q·η| ∈ (1, 2.4).
- **ε_vtx (vertex-pointing eff):** few-% inefficiency from the d0pair/z0sinθpair
  cuts, from simulation, validated on the UPC sample (cuts remove 1.4% / 0.5% of
  UPC pairs); MC d0 smeared to match data.
- **Net effect:** correction increases the yield by ≈2× at m_μμ = 8 GeV / p̄T = 4
  GeV, and by ≈30% (peripheral) to ≈50% (0–10% central) at higher m_μμ, p̄T.

**HF-decay background (§4.2) [background-for-writing — their bkg = our signal].**
After vertex pointing, the surviving HF background is estimated by a
**template fit to the d0pair distribution**: signal template from STARlight+HIJING
MC; background template from **data** (pairs with A > 0.06 ∧ α > 0.012, i.e. NOT
back-to-back — effectively pure background). Poisson log-likelihood fit (MINUIT)
over d0pair < 0.3 mm gives a fit signal fraction, translated to the analysis
range d0pair < 0.1 mm via Eq. (2). Signal fraction rises with p̄T and toward
peripheral: ≈50% in 0–5% (p̄T>4 GeV), →1 for UPC; quenching of HF in central
collisions weakens the centrality dependence at high p̄T. (This d0-based
separation is **complementary** to our Δp/p purity, [[muon_source_template_fits]].)

**DY background (§4.4).** Powheg+Pythia8 (AZNLO tune, CTEQ6L1), 5 nuclear PDFs
(**nNNPDF2.0 nominal**); only ≈1.8% / 1.2% of preselected DY survive Fid-α /
Fid-k⊥. Normalized per centrality via N_DY = (L·σ_PbPb_had·Δcent)·(σ^Fid_DY,NN·⟨T_AA⟩)
(Eq. 4). DY largest in central collisions (T_AA scaling). ~30% nPDF spread.

**Centrality (§3.3, Table 1) — 5.02 TeV ⟨T_AA⟩ [mb⁻¹]:**
0–5%: 26.0 · 5–10%: 20.4 · 10–20%: 14.4 · 20–30%: 8.77 · 30–40%: 5.09 ·
40–50%: 2.75 · 50–60%: 1.35 · 60–70%: 0.601 · 70–80%: 0.239 · 80–90%: 0.0815.
FCal ΣE_T with standard ATLAS Glauber calibration; ≥1 neutron in each ZDC required
to suppress UPC. (Our analysis: 5.36 TeV T_AA — different numbers; roadmap §Q2.3.)

**Cross-section & unfolding (§5.1).** σ_cent = (1/L)(N_sig − N_DY); "normalized
yields" Y = σ_cent/σ_tot cancel normalization systematics. Differential
cross-sections obtained from **unfolded** α, k⊥ distributions (D'Agostini
iterative Bayesian, RooUnfold; Refs. [54,55]).

**Systematic uncertainties (§5.2).** Luminosity **1.5%**; muon WP (medium→tight)
**7%** on cross-sections / 2% on normalized yields; trigger eff within **3%**
(rising to **6%** in 0–10%); reco eff **≈2%** crossx / <0.5% yields; d0pair /
z0sinθpair vertex cuts **≈2%** (0–5%); background-shape parameterization ≈0.2%;
d0pair signal-template variation **≈1%**; "excess background" (unidentified, ~4%
of fiducial yield) one-sided ≈4% crossx / 2% yields. **Total ≈8%** on
cross-sections. (Yields: ~69490 pairs Fid-α, ~67789 Fid-k⊥.)

**Result (Conclusion).** Cross-sections vary weakly with centrality; α and k⊥
distributions broaden toward central with a depletion near zero; no |Δy| or
event-plane dependence → **no significant magnetic-field effect** on the muons.

## References worth future reading   (≤3)

1. **D'Agostini iterative Bayesian unfolding / RooUnfold** — Refs. [54]
   (D'Agostini, Nucl. Instrum. Meth. A 362 (1995) 487) + [55] (Adye, RooUnfold,
   arXiv:1105.1160). *SUPPORTIVE.* New info: the standard **unfolding method**
   used here (and across ATLAS HI). Serves: detector response / **unfolding**
   (roadmap step 13 / overview §4c) — currently **no KB doc covers unfolding**.

Not added (already covered or superseded):
- [41] ATLAS Run 2 muon reco eff (arXiv:2012.00578) → already
  [[ATLAS_Run2_muon_reconstruction]] (the per-muon eff & data-MC SF this paper builds on).
- [49] HF-muon R_AA (arXiv:2109.00411) → already [[run2_hf_muon_raa]];
  [48] HF muons 2.76 TeV (arXiv:1805.05220) → already in FUTURE_READ via the
  back-to-back paper (Δp/p purity); [50] open-HF review (arXiv:1903.07709) →
  [[open_hf_production]]; [40] dense-env tracking (arXiv:1704.07983) →
  [[atlas_inner_detector_tracking]].
- [7] the **previous** γγ→μμ acoplanarity measurement (PRL 121 (2018) 212301,
  arXiv:1806.08708) — this paper **repeats and supersedes it** with 4× luminosity;
  not worth a separate summary.

## Related KB docs
- [[run2_dimuon_backtoback_paper]] — the back-to-back Letter that **cites this as
  Ref. [60]** for its per-pair trigger×reco efficiency weight; this doc supplies
  the cited method definitions.
- [[run2_dimuon_note]] — our reference internal note, where **γγ→μμ is a
  background removed** from OS pairs (the α/A region this paper selects).
- [[ATLAS_Run2_muon_reconstruction]] — Ref. [41]: single-muon reco/ID efficiency
  and data-MC scale factors this paper's reco efficiency builds on.
- [[atlas_run2_muon_trigger]] — single-muon trigger-efficiency (tag-and-probe,
  dimuon factorization) methodology underlying the per-pair trigger combination here.
- [[muon_source_template_fits]] — our Δp/p purity framework; this paper's
  d0pair HF/signal template fit is a complementary d0-based separation.
- [[run2_hf_muon_raa]] — Ref. [49]: HF-muon R_AA (their dominant background = our
  signal); HF quenching context.
- [[open_hf_production]] — Ref. [50] hub: HF energy loss / quenching that
  suppresses their central-collision HF background.
- [[atlas_inner_detector_tracking]] — Ref. [40]: dense-environment tracking and
  close-track ambiguity that **justifies OUR pair ΔR axis** (and which this
  paper's back-to-back pairs avoid, enabling its ε₁·ε₂ factorization).
</content>
</invoke>
