# Run 2 Dimuon Back-to-Back Δφ Correlation — Published Letter (PRL)

**Source:** ATLAS Collaboration, *Azimuthal angle correlations of muons produced
via heavy-flavor decays in 5.02 TeV Pb+Pb and pp collisions with the ATLAS
detector* — Phys. Rev. Lett. **132** (2024) 202301; CERN-EP-2023-185.
**arXiv / DOI:** 2308.16652 [nucl-ex] / 10.1103/PhysRevLett.132.202301
**PDF:** `./2308.16652.pdf` (PRIMARY)
**Classification:** PRIMARY
**Added:** 2026-06-14

> This is the **published Letter** of the **same analysis** summarized in
> [[run2_dimuon_note.md]] (internal note ATL-COM-PHYS-2021-1094). The note has the
> implementation detail (exact cuts, fit machinery, systematics tables); **this
> doc captures the publication-level framing** — physics questions, comparisons,
> implications — and the few numbers/methods stated cleanly in the Letter. Use
> both together; do not duplicate. Cross-linked from the note's header.

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Per-pair weight = inverse product of trigger × reco efficiencies, ε(pT,η) of the two muons; PbPb uses two dimuon triggers, pp uses the 2mu4 trigger only | Trigger-eff (roadmap step 8 PbPb single-mu4 / step 10 pp24 2mu4; IntNote §8). **We transfer the per-pair trigger-eff method but ADD a ΔR correlation factor** (overview §4b) | [method-we-use] |
| Combinatorial/underlying-event background built from same-sign (SS) pairs; cos(2Δφ) elliptic modulation of the combinatoric pedestal | Combinatorial subtraction via SS pairs (overview §4a; IntNote yield section). The UE/combinatorial estimation transfers directly | [method-we-use] |
| Tight muon WP keeps π/K-decay fraction < 5%; medium WP used as systematic | Muon WP choice (roadmap Q2.8), fake-muon purity context (step 16) | [method-we-use] |
| Δp/p-style π/K purity assessed via procedure of Ref. [19] | Δp/p significance template-fit purity framework (roadmap step 16, BLOCKED) | [method-we-use] |
| OS-only pair selections: Υ [9.2,10.4] & Z [70,110] GeV mass vetoes; γγ→μμ removal; |Δη|>0.8 quarkonia/jet-peak suppression | Pair-selection definitions (IntNote §7); acceptance-correction concept | [method-we-use] |
| DY background subtraction via Powheg+Pythia8 (nNNPDF2.0; syst nCTEQ15/NNPDF3.0) | DY treatment if/when needed (OS pairs); IntNote backgrounds | [background-for-writing] |
| Big physics question: angular deflection / collisional-vs-radiative scattering of b in QGP; novel first dimuon-from-HF angular-correlation measurement; model-independent constraint | Intro / motivation phrasing for IntNote & thesis; HF-energy-loss narrative | [background-for-writing] |
| Comparison to prior HF probes (HF-hadron & decay-lepton R_AA and v_n; b-jet quenching) | Intro citations / where our R_AA sits in the field | [background-for-writing] |

## Scope & condition-difference warnings

Run 2, √s_NN = **5.02 TeV**: Pb+Pb **1.94 nb⁻¹** (2015+2018), pp **0.26 fb⁻¹**
(2017) (Letter p.2). Our analysis is **Run 3, 5.36 TeV** (pp24 + PbPb 23/24/25).
ACKNOWLEDGE the energy/conditions/detector difference; DO NOT assume its size or
direction.

**Signal-topology difference (do not quantify physics implications).** This Letter
measures the **away-side (Δφ ≈ π) back-to-back** correlation between muons from
**two different** heavy quarks created in one hard scatter (momentum conservation
→ back-to-back). **Our** signal is the opposite: a **single-b** pair where both
muons come from **one** b-decay chain (b→μ and b→c→μ), giving **nearby / small
opening angle**, low mass. The Letter's away-side observable (Γ, σ of the Δφ peak)
is NOT our observable; our observables are differential cross-sections and R_AA
(overview §3). Whether a Δφ-width measurement is also in our scope is open
(roadmap §Q2.9).

**Reco-efficiency difference.** This Letter weights by a **product of two
single-muon** reco efficiencies ε_reco(pT,η)·ε_reco(pT,η). **Our** reco efficiency
is a single **pair** efficiency ε_reco(pair pT, pair η, ΔR) — because our two
muons are close-by and reconstruction-correlated, it is **not** ε₁·ε₂ (overview
§4b; [[../analysis/run2_dimuon_note.md]] differs likewise). The trigger-eff
method transfers; the reco-eff is genuinely a different object.

## Content summary

**Physics motivation (Letter pp.1–2).** HF quarks are produced early in hard
scatters and probe the full QGP evolution; their interaction with the medium
causes energy loss (suppression of HF hadrons / decay leptons / b-jets) and flow
(HF v_n). A complementary probe is the **angular deflection** of the quarks:
quark–medium scattering should deflect quarks (Brownian diffusion at p≲m_Q; at
p≫m_Q collisional broadening, damped by radiative gluon emission). Theory (Ref.
[44], Nahrgang et al., PRC 90 (2014) 024907) shows **b-quark angular correlations
are sensitive to the relative weight of collisional vs radiative scattering**.
Direct b-hadron-pair detection is hard; **muon pairs from simultaneous
semileptonic B decays are feasible — no such measurement existed before this
Letter** (Letter p.2).

**Charge handle (Letter p.2).** Both charm and bottom give single muons, but for
**both muons with pT > 4 GeV** the charm dimuon contribution is **kinematically
suppressed by ~an order of magnitude** vs bottom (Powheg+Pythia8, Ref. [51]
supplemental). cc̄ → almost exclusively OS muon pairs; bb̄ → both OS and SS (b→c
sign reversal + neutral-B mixing, ~half SS when a parent is a neutral B). Hence
**SS pairs ≈ clean bb̄ probe**; OS−SS comparison estimates the cc̄ contribution.
In pp, bottom pairs contribute **~90%** of the OS yield (Letter summary).

**Selections (Letter pp.2–3).** Single muon: **tight** WP, pT > 4 GeV, |η| < 2.4.
Pair: average **p̄T = (pT1+pT2)/2 > 5 GeV** (removes low-efficiency pairs, improves
precision); **|Δη| > 0.8** (suppresses quarkonia & in-jet HQ pairs, removes light
mesons & J/ψ, but admits Υ). OS-only extra cuts: mass vetoes **[9.2,10.4] GeV**
(Υ) and **[70,110] GeV** (Z); **γγ→μμ removal** keep if |Δφ−π|/π > 0.01 OR
|pT⁺−pT⁻|/(pT⁺+pT⁻) > 0.08. Multi-muon events (~2%): all pairs used. π/K-decay
fraction < 5% (tight WP, via Ref. [19] procedure).

**Triggers (Letter p.2).** Two dimuon triggers in Pb+Pb: (1) single μ pT>4 at L1
+ two μ pT>4 at HLT; (2) two μ pT>4 at both L1 and HLT. pp used only trigger (2).
(Further detail in Ref. [60].)

**Efficiency correction (Letter p.3).** Each pair weighted by inverse product of
**reconstruction and trigger** efficiencies, computed vs pT and η of the two
muons (method in Ref. [60]). Average weight ≈ **2.3 (pp), 2.4 (Pb+Pb)**. OS pairs
get an extra **acceptance correction** ε_opp(Δφ) for the OS-only mass/γγ cuts,
obtained by applying those cuts to SS pairs (default) cross-checked with
mixed-event OS pairs (single-muon-triggered sample); difference < 0.5% → syst.

**Background.** Combinatorial = pedestal C_comb (geometrically enhanced in Pb+Pb;
modulated by cos(2Δφ) from QGP flow). DY subtracted via Powheg+Pythia8 nNNPDF2.0
(syst nCTEQ15, NNPDF3.0).

**Observable & fit (Letter pp.3–4, Eqs. 1–2).** Correlation function
C(Δφ) ≡ (1/N_tot) dN/dΔφ, 32 bins over [−π/2, 3π/2], symmetrized about 0.
Fit: `C_fit(Δφ) = C_comb[1 + 2 v_2,2^eff cos(2Δφ)] + C_corr(Δφ)`, with
`C_corr = C_max Γ² / ((Δφ−π)² + Γ²)` (Cauchy-Lorentz, folded at −π/2 and 3π/2).
Width observables: **Γ** = HWHM, and **σ** = std dev of C_corr around Δφ=π with
the pedestal C_corr(0) subtracted. Cauchy-Lorentz describes the peak well;
Gaussian/gen-Gaussian/von-Mises fail. ∫C_corr(OS) ≈ 2× ∫C_corr(SS) (matches
Pythia8). v_2,2^eff ~ O(0.01) in Pb+Pb, equal for SS/OS (→ flow of uncorrelated
HF muons). Stat errors via Gaussian re-sampling + refit.

**Systematics (Letter p.5).** Sources: muon WP (tight↔medium), trig+reco eff ±1σ,
OS mass/γγ acceptance (same-event SS vs mixed-event), combinatoric
parameterization, width-extraction method. pT/η/φ resolution negligible.
**Leading systematic = combinatoric parameterization** (add cos3Δφ for Pb+Pb; for
pp constrain v_2,2^eff→0): **~1% (2.5%)** on σ and **~2% (5.5%)** on Γ in Pb+Pb
(pp).

**Results (Letter pp.6–7, Figs. 1–2).** Clear away-side Δφ≈π peak in all samples.
Widths (σ, Γ) show **no significant centrality variation over 10–80%** and are
**consistent with pp**; a **significant decrease (narrowing) in 0–10% central**
(~2σ for all-pairs). SS and OS widths similar (bb̄ dominance). The QGP
extra-broadening σ_int² = σ²_PbPb − σ²_pp is consistent with zero outside 0–10%;
**90% CL upper limit σ_int ≲ 0.2** (all-pairs). Theory expectation (Ref. [44]) of
broadening increasing toward central is **not supported**. The ⟨p̄T⟩ is
centrality-independent except ~3% higher in 0–10% (pp ⟨p̄T⟩ = 7.36 ± 0.10 GeV over
5–20 GeV, Ref. [73]); Powheg+Pythia8 b(c)-hadron-to-muon angular rms width =
**0.158 (0.058)**, ≪ muon-pair width → decay smearing small. **Conclusion:** b
quarks in QGP are strongly quenched (energy loss) but suffer little angular
deflection — a **model-independent constraint** on stochastic b-quark deflection.

**Supplemental (Letter pp.10, Figs. 3–4).** Powheg+Pythia8 bb̄→μμ + cc̄→μμ
(nCTEQ15+SIH nPDF), same fiducial cuts, reproduces the data widths; OS sample
~85–96% bb̄ depending on charge combo.

## References worth future reading   (≤3)

1. **ATLAS, γγ→μμ in non-UPC Pb+Pb at 5.02 TeV** — PRC 107 (2023) 054907,
   arXiv:2206.12594 (Ref. [60]). *PRIMARY.* New info: the **actual definition of
   the per-pair trigger & reconstruction efficiency method** and the γγ→μμ removal
   that this Letter only cites. Serves: trigger-eff (steps 8/10) and reco-eff
   (step 12) method provenance, pair selection.
2. *(Already in the KB — see [[run2_hf_muon_raa]]; Ref. [21] = arXiv:2109.00411,
   the single-muon HF R_AA. Not a future-read.)*
3. **ATLAS, suppression & azimuthal anisotropy of HF-decay muons in Pb+Pb 2.76
   TeV** — PRC 98 (2018) 044905, arXiv:1805.05220 (Ref. [19]). *SUPPORTIVE.* New
   info: the **π/K fake-muon purity procedure** (Δp/p) that both this Letter and
   our analysis adopt. Serves: Δp/p purity framework (step 16).

(Not added: Ref. [44] Nahrgang theory — relevant only to the away-side
Δφ-broadening interpretation, which is not our nearby-pair observable.)

## Related KB docs
- [[run2_dimuon_note.md]] — same analysis, internal note: full implementation
  detail (exact cuts, fit code, systematics tables). This doc adds the published
  framing/comparisons it lacks.
- [[../physics/detector/ATLAS_Run2_muon_reconstruction.md]] — muon reco/ID
  performance & tag-and-probe efficiency methodology underlying the per-muon
  efficiencies used here.
- [[../../../Analysis/docs/analysis_overview.md]] — our Run 3 single-b signal,
  pair reco-eff (pair pT, pair η, ΔR), and trigger-eff (with ΔR factor) methods
  that adapt this Letter's machinery.
</content>
</invoke>
