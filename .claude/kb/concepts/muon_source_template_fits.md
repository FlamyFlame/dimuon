# Concept: muon-source template fits (Δp/p momentum-imbalance & d0 impact-parameter)

**Type:** cross-paper **concept hub** (not a source summary). Synthesizes a method
that recurs across three KB source docs; each source doc keeps only its
specificity and links here. Per `KB_BUILDING_GUIDE.md` §8 (dedup) and §5
(attribute; preserve genuine per-paper differences — these fits are used
*differently* in different papers, see below).
**Added:** 2026-06-15

> **Two distinct fits, often used together.** (1) the **Δp/p momentum-imbalance**
> fit separates **real muons from hadronic/fake (π/K) muons**; (2) the **d0
> impact-parameter** fit separates **charm from bottom**. They are statistically
> independent (signal ρ⊥d0), so they factorize into a two-step 1-D fit.

## Relevance to this analysis
- **Δp/p (significance) template fit = OUR fake-muon purity framework**
  (`analysis_overview` §6, roadmap step/`task_07`, IntNote §16 "Background /
  fake-muon purity"). `analysis_overview` states the method is "adapted from the
  Run 2 note." `[method-we-use]`.
- **d0 charm/bottom fit = contrast, not our method.** We separate b from c by the
  **single-b dimuon kinematics** (low mass, small ΔR) + SS/OS, not by a per-muon
  d0 fit. Keep d0 as `[background-for-writing]` (Intro/method comparison).

## (1) Δp/p momentum-imbalance fit — real vs fake muon
- **Variable:** ρ = (p_ID − p_MS)/p_ID, with p_MS the muon-spectrometer momentum
  energy-loss-corrected back to the ID. **Real** muons (HF + prompt) give a
  symmetric peak at ρ ≈ 0; **hadronic/fake** muons (π/K decay-in-flight,
  punch-through, random ID+MS combinations) are broad and **shifted to positive ρ**
  (their MS momentum is mismodeled).
- **Templates:** muon-filtered / non-diffractive Pythia8 QCD; muon momentum
  **shifted & smeared** to match data, tuned on the J/ψ→μμ invariant-mass response
  (vs centrality in Pb+Pb). Pb+Pb ρ templates use **no overlay** (overlay ID
  resolution worse than data).
- **Two ways it is used in the KB sources (the key per-paper difference):**
  - **Yield extraction** — [[run2_hf_muon_raa]], [[run2_hf_muon_vn]]: a
    multi-component template fit (signal HF, hadronic, fake, prompt-fixed) returns
    the **HF-muon yield** separated from hadronic+fake background. This is a
    correction that enters the measured yield.
  - **Purity demonstration** — [[run2_dimuon_note]]: ρ is converted to a
    **significance** (ρ−μ(pT,qη))/σ(pT,qη) using true-muon μ,σ from HIJING-overlay
    MC (making it pT/η/centrality-independent); the **pair significance** =
    quadrature sum of the two muons; the data distribution is fit with **three**
    templates (Sig–Sig, Sig–Bkg, Bkg–Bkg). Result: **>98% purity → NO cut applied**
    (a purity demonstration, not a yield correction). **This is the version OUR
    analysis adapts.**

## (2) d0 impact-parameter fit — charm vs bottom
- **Variable:** muon track transverse impact parameter d0 (wrt beam spot / primary
  vertex). Different **B vs D lifetimes** (B≈1.5, D⁺≈1.0, D⁰≈0.4 ×10⁻¹² s) give
  different d0 spreads (bottom broader). Bottom includes the b→c→μ cascade.
- **Fit:** two-component over |d0| < 0.5 mm; single free parameter = **bottom
  fraction f_b = κ_b/(κ_c+κ_b)**, applied *within* the HF yield fixed by the ρ fit.
- **Templates:** muon-filtered Pythia8 reweighted to **FONLL** c/b spectra (+
  baryon/meson, D⁺/D⁰ corrections); in Pb+Pb additionally reweighted to measured
  (ALICE/CMS) modified hadron spectra and built with **MB Pb+Pb overlay** (for
  vertex resolution). Used in [[run2_hf_muon_raa]], [[run2_hf_muon_vn]]; **not**
  used in [[run2_dimuon_note]] (kinematic c/b separation there).

## Provenance / which source to read for detail
- Origin of the d0 c/b method: ATLAS PRL 124 (2020) 082301 (arXiv:1909.01650, pp)
  — in FUTURE_READ.
- Yield-fit detail + systematics: [[run2_hf_muon_raa]] (note §4.4), [[run2_hf_muon_vn]] (note §6).
- Significance/purity (our version) detail: [[run2_dimuon_note]] (Δp/p section).

## Related KB docs
- [[run2_hf_muon_raa]], [[run2_hf_muon_vn]], [[run2_dimuon_note]] — the three sources using these fits.
- [[atlas_gammagamma_mumu_pbpb]] — uses a **d0pair** template fit (MC signal vs data-sideband) for HF separation — a complementary d0-based alternative to the Δp/p purity here.
- [[yun_tian_thesis]] — Run-1 (2.76 TeV) origin of the Δp/p "Eloss" momentum-balance template fit, with full method depth.
- [[../physics/detector/ATLAS_Run2_muon_reconstruction]] — d0/z0 resolution and the ID/MS momentum measurements these fits exploit.
- `Analysis/docs/analysis_overview.md` §6 — our fake-muon purity program (ground truth).
</content>
