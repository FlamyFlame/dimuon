# Run 2 Dimuon Internal Note — Comprehensive Summary

**Reference:** `~/workarea/dimuon_codes/IntNotesRun2DimuonReference/`
(ATL-COM-PHYS-2021-1094, ANA-HION-2022-19). Compiled PDF:
`IntNotesRun2DimuonReference/ATL-COM-PHYS-2021-1094.pdf` (51 pp).
LaTeX sources in `IntNotesRun2DimuonReference/tex/*.tex` (read those, not
the PDF — no poppler on SDCC nodes).

**Title:** *Measurements of the azimuthal correlations of muon pairs from
decay of heavy-flavor quarks in Pb+Pb and pp collisions* (Cole, Krauth,
Mohapatra — Nevis/Columbia).

This is the completed Run 2 (5.02 TeV) analysis that the current Run 3
(5.36 TeV) analysis follows. Use this file for structure/methodology
context; consult the tex sources for exact wording, equations, figures.

**Published as PRL 132 (2024) 202301 (arXiv:2308.16652).** For the
publication-level framing (physics questions, field comparisons,
implications) see the companion doc [[run2_dimuon_backtoback_paper.md]];
this note is the implementation detail behind it.

## Physics goal

Measure Δφ (azimuthal separation) correlations of HF-decay muon pairs in
PbPb and pp. Extract away-side (Δφ=π) correlated **widths** (and
optionally yields — note had an `IncludeYields` LaTeX toggle, default
OFF) vs centrality. QGP interactions of HF quarks are predicted
(Nahrgang et al.) to broaden the away-side peak; collisional vs radiative
energy loss give different broadening. **Headline Run 2 result:** widths
are centrality-independent and consistent with pp (slight *narrowing* in
0–10% central) → HF quarks are strongly quenched but suffer little
angular deflection. 90% C.L. upper limits set on QGP broadening
σ_int = √(σ²_PbPb − σ²_pp).

SS pairs ≈ pure bb̄ (B-oscillation + b→c cascade); OS pairs from bb̄ and
cc̄. Pythia studies (note's Muon Sources appendix): with 2μ pT>4 GeV
|η|<2.5 requirement, bb̄ away-side yield ~10× cc̄; OS:SS ≈ 2:1 for bb̄.

## Publication-level framing (added from PRL 132 (2024) 202301)

Phrasing/comparisons the published Letter adds beyond the note's terse intro
(use for IntNote/thesis motivation; see [[run2_dimuon_backtoback_paper.md]]):

- **Novelty:** first measurement of HF angular correlations via **pairs of
  muons from simultaneous semileptonic B decays** — direct b-hadron-pair
  detection is hard, so the dimuon channel is the feasible probe.
- **Theory framing:** quark–medium scattering deflects HF quarks — Brownian
  diffusion at p≲m_Q, collisional broadening at p≫m_Q damped by radiative
  gluon emission; b-quark angular correlations are sensitive to the
  **collisional-vs-radiative** balance (Nahrgang et al., PRC 90 (2014)
  024907). Away-side Δφ width is the observable for this.
- **Charm handle (quantified):** for both muons pT>4 GeV, charm dimuons are
  **~an order of magnitude** below bottom (Powheg+Pythia8); bottom pairs are
  **~90%** of the OS yield in pp; SS pairs ≈ clean bb̄.
- **Implication:** widths consistent pp↔PbPb and centrality-independent
  (10–80%), narrowing only in 0–10% → b quarks strongly quenched but little
  angular deflection; **model-independent 90% CL limit σ_int ≲ 0.2** on
  extra QGP-induced broadening. Theory expectation of central broadening is
  not supported.
- **Field context:** complements HF-hadron/decay-lepton **R_AA** and **v_n**
  and b-jet quenching measurements (cite ATLAS HF-muon R_AA PLB 829 (2022)
  137077 and v_n PLB 807 (2020) 135595).

## Datasets (Run 2)

- **PbPb 5.02 TeV:** 2015 (0.50 nb⁻¹, 33 runs) + 2018 (1.44 nb⁻¹, 39
  runs) = 1.94 nb⁻¹. GRLs listed in `Data_and_selections.tex`. 2018 GRL
  excludes runs 367318–367384 (toroid issue). HardProbes stream AODs.
- **pp 5.02 TeV:** 2017 reference run, 0.26 fb⁻¹, μ⟨pileup⟩ 0.4–4,
  `physics_Main` periodM.
- **Centrality:** FCal ΣET, HI-group standard bins, 0–80% (80% ⇔
  ET_FCal = 63.7 GeV). ⟨T_AA⟩ from HI group Glauber (table in
  `DphiFits.tex`).

### Triggers
- PbPb signal: `HLT_mu4_mu4noL1` (primary; prescaled in a few runs) OR
  `HLT_2mu4` (unprescaled; not a subset because mu4 vs mu4noL1 HLT hypos
  differ). Combination rule: if pair matches mu4_mu4noL1 → use its
  efficiency; else if matches 2mu4 AND mu4_mu4noL1 was prescaled out in
  that event (xAOD flag) → use 2mu4 efficiency; else drop pair (counted
  as mu4_mu4noL1 inefficiency).
- PbPb trig-eff support: MinBias (TE600, TE50_VTE600, ZDC sptrk VTE50
  triggers) + single-mu HLT_mu4/mu6/mu8 (HardProbes).
- pp signal: `HLT_2mu4`, unprescaled.

### MC samples
- **PbPb reco eff + signal Δp/p templates:** STARlight γγ→μμ + HIJING
  *overlay* MC (made for ANA-HION-2020-10; JIRAs ATLHI-304 2015-cond,
  ATLHI-313 2018-cond). 2015/2018 conditions combined.
- **PbPb background Δp/p templates:** Pythia8 jetjet JZ1–4 HIJING-overlay
  (mc16_valid 420011–420014).
- **pp:** reco/trig eff taken from MCP recommendations / HION-2019-58
  (not re-derived). Pythia8 dijet JZ0–5 mc16_5TeV samples for Δp/p
  templates only.
- **DY correction:** official Powheg DY sample (ATLHI-331).

## Muon & pair selections

Single muon: combined; **Tight** WP default (Medium for syst/cross-check);
pT > 4 GeV; |η| < 2.4; |d0| < 5 mm.
Pair: trigger-matched (TriggerMatchingTool, dR<0.1); **p̄T ≡ (pT1+pT2)/2
> 5 GeV** (S/B); all pairs used in multi-pair events (~2%).
Pair-level cuts (PbPb *and* pp, for identical acceptance):
- **γγ→μμ removal (OS only):** reject Aco<0.01 && A_μμ<0.08
  (Aco = 1−|Δφ|/π; A = |pT⁺−pT⁻|/(pT⁺+pT⁻)). <1% effect in pp.
- **Mass vetoes (OS only):** remove M_μμ 0.6–1.05, 2.9–3.3, 3.6–3.8,
  9.2–10.4, 70–110 GeV (ρ/ω/φ, J/ψ, ψ′, Υ, Z).
- **|Δη| > 0.8:** removes near-side (0,0) jet peak; also kills all
  resonances below M~3.4 GeV.

## Corrections

Per-pair weight w⁻¹ = ε_trig(p_a,p_b) · ε_reco(p_a) · ε_reco(p_b).

1. **Trigger eff (PbPb):** factorized from single-muon efficiencies.
   ε(2mu4) = ε_mu4(a)·ε_mu4(b). ε(mu4_mu4noL1) = inclusion–exclusion
   with mu4 and mu4noL1 legs incl. overlap term (eq. in
   `Corrections.tex`). ε_mu4 from MinBias tag-free matching (offline μ
   matched to HLT μ dR<0.1) vs (pT, q·η); fitted Fermi + linear-tail,
   data-point interpolation used below pT=8 GeV. ε_mu4noL1 and
   mu4∩mu4noL1 from HardProbes single-mu stream: tag = first muon
   matching mu4/mu6/mu8, probe fraction vs (pT, q·η) of second muon.
   Cross-check: mu4 eff from HP stream = MinBias result, validating the
   factorization. Centrality dependence applied as a multiplicative
   factor CentDep(pT, centrality) (not enough stats for 3-D).
2. **Reco eff (PbPb):** from STARlight+HIJING-overlay MC; ε = fraction of
   generated muons with reco match (match-prob>0.5) passing analysis
   cuts, vs (pT, q·η, centrality 10–20% bins); data-point interpolation
   (no fits); MCP scale factors on top.
3. **pp:** reco + trig eff from prior publication (HION-2019-58) / MCP.
4. **M_μμ / γγ-cut Δφ-acceptance correction (OS only):** the OS-only mass
   and Aco/Asym vetoes sculpt Δφ. Efficiency vs Δφ estimated (a) default:
   applying cuts to SS pairs; (b) alternate: event-mixing with HLT_mu4
   single-muon sample, centrality-matched. OS pairs reweighted by 1/eff.
   Difference (a)-(b) → systematic.
5. **Drell-Yan subtraction (OS only):** dN_DY/dΔφ = L·σ_had^PbPb·Δcent·
   ⟨T_AA⟩·dσ_DY,NN/dΔφ from Powheg (ATLHI-331); nNNPDF2.0 default,
   nCTEQ15/NNPDF3.0 variants (check only).

## Background (fake-muon) estimate — Δp/p template fits

> **Shared mechanism (ρ=Δp/p real-vs-fake) lives in the concept hub
> [[../concepts/muon_source_template_fits]] — this is the OUR-adapted purity
> variant.** Below = *this analysis's specifics* (significance transform, pair
> significance, 3-template purity demonstration, >98% → no cut).

Momentum imbalance Δp/p = (p_ID − p_MS,extrapolated)/p_ID. HF muons:
Gaussian around 0; π/K in-flight decays & calo secondaries: broad,
shifted positive. Converted to a **significance** (Δp/p − μ(pT,qη))/
σ(pT,qη) with μ,σ for true muons from HIJING-overlay MC → cut/template
efficiency independent of pT, η, centrality. Pair significance =
quadrature sum of the two muons. Data pair-significance distribution fit
with 3 MC templates: Sig–Sig, Sig–Bkg, Bkg–Bkg. **Result: >98% purity →
no Δp/p cut applied in the analysis** (the study is a purity
demonstration, not a correction).

## Signal extraction — Δφ fits

dN/dΔφ = N₀(1 + v_2,2 cos 2Δφ) + Y^corr(Δφ),
Y^corr = N^corr/((Δφ−π)² + Γ²) − N^pedestal (Lorentzian + ZYAM pedestal
forcing Y^corr(0)=0; folded into (−π/2, 3π/2)).
- N₀: flat UE; v_2,2 term: elliptic-flow modulation of UE.
- **Width observable = RMS of Y^corr** (remapped so peak is centered),
  NOT Γ. Γ also reported.
- Statistical errors via toy resampling+refit (StatErrs appendix).
- Fits per centrality bin, per charge combo (SS / OS / All), Tight &
  Medium, PbPb and pp.

## Systematics (PbPb; pp analogous in SystematicsPP appendix)

| Source | Method | Width impact (SS/OS/All %) |
|---|---|---|
| Muon WP | Tight↔Medium | 2 / 0.5 / 0.5 |
| Reco+trig eff | vary ±1σ (MCP / HION-2020-10) | 0.1 |
| Δη cut | 0.8→0.9 | 0.5 / ~0 / 0.25 |
| Fit model | fold-subtract alternative vs nominal | ~0 / 2.0 / 0.75 |
| Higher-order vₙ | add cos3Δφ | 1.5 / 1 / 1 |
| M_μμ correction | same-event SS vs event-mixed eff | ~0 |
| Drell-Yan PDF | 3 PDFs (check only) | 0.5 / 0.3 (OS/All) |
| Trigger choice | 2mu4-only (check only) | — |

Smoothed by fitting %-variation vs centrality; total = quadrature sum.
Lumi (±1.5% PbPb, ±1.6% pp) & T_AA only matter for yields (toggle off).

## Note structure (template for Run 3 note)

Intro → Datasets (data, triggers, GRL, MC) → Efficiency Corrections
(trig, reco, MC closure) → Analysis (selections, observables, M_μμ/γγ
acceptance correction, DY correction, Δφ fits, widths) → Momentum
Imbalance & Template Fitting (purity) → Systematics → Results/Summary.
Appendices: fit params, pp systematics, Medium-muon plots, stat-error
toys, Muon Sources (Pythia origin study), reco-eff centrality dep.,
mean-pT spectra, multi-pT-interval fits, linear-fit slopes, Powheg
comparisons. (Unused/commented: separate Templates, DrellYan,
DphiFitsDeta0p9 sections.)

## Key differences to keep in mind for the Run 3 analysis

- Run 3: 5.36 TeV PbPb 2023/24/25 + pp24 reference; Run 2: 5.02 TeV.
- Run 3 PbPb uses single-mu4 trigger (trigger_mode=1) as primary, pp24
  uses 2mu4 — different trigger-efficiency machinery (see
  `Analysis/docs/tracking/mu4_trig_effcy_implementation.md`).
- Run 3 reco eff from Pythia hQCD fullsim + HIJING overlay (not
  STARlight overlay); pp eff from Pythia fullsim pp24.
- Run 2 measured widths only (yields toggled off); Run 3 scope includes
  cross-sections / R_AA (crossx pipelines).

## Related KB docs   (knowledge graph)
- [[run2_dimuon_backtoback_paper]] — the published Letter of THIS analysis (publication framing).
- [[../concepts/muon_source_template_fits]] — the Δp/p significance / pair-significance purity fit (our adapted version) lives here as a concept hub.
- [[run2_hf_muon_raa]], [[run2_hf_muon_vn]] — same-group HF-muon analyses sharing the Δp/p method + FCal centrality.
- [[physics/heavy_ion/open_hf_production]] — field background (energy loss, R_AA, dead-cone).
- [[physics/background/gluon_splitting_flavour_excitation]] — the g→QQ̄ correlated background and the charm-vs-bottom dimuon handle.
