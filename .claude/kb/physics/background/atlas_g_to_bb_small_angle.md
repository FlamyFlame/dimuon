# ATLAS g→bb̄ at small opening angle (pp 13 TeV)

**Source:** ATLAS Collaboration, *Properties of g→bb̄ at small opening angles in pp collisions with the ATLAS detector at √s = 13 TeV* — Phys. Rev. D 99 (2019) 052004; CERN-EP-2018-323
**arXiv / DOI:** arXiv:1812.09283 ; doi:10.1103/PhysRevD.99.052004
**PDF:** `./1812.09283.pdf` (PRIMARY)
**Classification:** PRIMARY
**Added:** 2026-06-15

> **Why this doc exists.** This is the direct citeable anchor flagged by the
> background-concept hub [[gluon_splitting_flavour_excitation]]: a dedicated ATLAS
> measurement of **g→bb̄ at SMALL OPENING ANGLE** — exactly our nearby-pair regime.
> It probes the very fragmentation (Pythia/Herwig/Sherpa parton shower of g→bb̄)
> whose modeling our Pythia g→QQ̄ background template inherits, and it quantifies
> the **large theory uncertainty** in the collinear region. This doc keeps the
> measurement specifics; the *mechanism* (FSR g→QQ̄ vs ISR flavour excitation, and
> signal-vs-background distinction) lives in the hub — do not duplicate it here.

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| **Direct differential measurement of g→bb̄ in the small-opening-angle / collinear regime** (ΔR(b,b) ∈ [0.2,1.0]), the kinematic regime of our correlated single-b-mimicking background (two muons from two close b's of one gluon split) | Intro / background framing; `analysis_overview` §6; IntNote background section | [background-for-writing] |
| **Data-vs-parton-shower-MC comparison (Pythia 8 / Sherpa 2.1):** Pythia mismodels the small-angle g→bb̄ fragmentation shapes (ΔR, log(mbb/pT), z(pT), Δθ); Sherpa generally closer but both disagree at low mass / low z / all Δθ | Validation/caution for our **Pythia g→QQ̄ background template** — the same shower modeling we rely on; motivates a template-modeling systematic — roadmap step 16 (background/purity), `task_07` | [method-we-use] (template validation) |
| **Large theoretical-modeling uncertainty in the collinear g→bb̄ region** (2–13% on these jet observables; dominant systematic alongside JES & unfolding) | Why our g→bb̄ background template carries a large modeling uncertainty; sets expectation that the collinear region is poorly constrained | [background-for-writing] |
| **Flavour-fraction fit to per-track signed impact-parameter significance s_d0** to separate BB / B / L+C jets (data-driven double-b purity) | Conceptual cousin of our d0/Δp/p purity template fits ([[muon_source_template_fits]]) — but on jet tracks, NOT muons; method contrast only | [background-for-writing] |

## Scope & condition-difference warnings
- **System:** pp, **√s = 13 TeV**, 2016 data, **33 fb⁻¹**, ⟨pileup⟩ ≈ 25 (p.4).
- **Probe is JET-level b-tagging, not muon-tagged dimuons.** g→bb̄ is tagged via
  two **R=0.2 track-jets** (b-quark proxies) inside one **R=1.0 trimmed calo jet**
  (gluon proxy) with double-b-tagging (MV2c10). Our analysis tags the same physics
  via **two muons** from semileptonic b decays. The OBSERVABLES (jet ΔR, jet z,
  jet mass) are not our observables (muon-pair m, pT, ΔR).
- **Not our beam energy / system:** Run 2 pp 13 TeV inclusive high-pT jets; our
  analysis is Run 3 **5.36 TeV pp and Pb+Pb**. ACKNOWLEDGE the difference; **DO NOT
  estimate the size or direction** of any transfer of rates/shapes to our
  conditions (GUIDE §5). What transfers conceptually is the *qualitative lesson*:
  parton-shower g→bb̄ modeling in the collinear region is imperfect and carries
  large theory uncertainty.
- **High-pT regime:** gluon (large-R jet) pT > 450 GeV — much harder than our
  signal-region pair pT > 8 GeV. The collinear small-ΔR enhancement here is driven
  by the high boost; do not assume the same shapes at our scale.

## Content summary

**Goal & topology (Intro, §4–5).** Measure differential properties of high-pT
g→bb̄ fragmentation at small opening angle — the dominant background to boosted
H→bb̄ searches, and a clean QCD probe (the b mass screens the soft-collinear
splitting singularity). Gluon proxy = anti-kt **R=1.0** calorimeter jet, trimmed
(reclustered R_sub=0.2, remove subjets with momentum fraction f_cut < 0.05).
b-quark proxies = anti-kt **R=0.2 track-jets** (tracks pT > 500 MeV, ≥1 pixel + ≥6
silicon hits, ≥2 tracks), ghost-associated to the large-R jet.

**Event selection (§4.2).** Single-jet trigger (fully efficient for offline
pT > 450 GeV, |η| < 2). Require leading calo jet pT > 450 GeV, ≥2 associated
track-jets with pT > 10 GeV, |η| < 2.5; leading two track-jets used. Leading
track-jet **b-tagged with MV2c10 at the 60% WP** (working point: ε_b = 60%,
ε_c = 15%, light-jet rejection 1/ε_light = 480, in tt̄ for jet pT > 10 GeV,
|η| < 2.5). Only ONE of the two is required b-tagged at detector level (requiring
both increases purity but degrades the background-fit precision). Particle level:
large-R jet pT > 450 GeV, both associated particle-track-jets tagged as b-jets
(b-hadron pT > 5 GeV ghost-matched).

**Observables (§5), unfolded to particle level (Fig. 6):**
- **ΔR(b,b)** between the two track-jets — opening angle; by construction
  0.2 ≤ ΔR ≤ 1.0, with a peak ~0.3 from the R=0.2 track-jet radius.
- **z(pT) = pT,2 / (pT,1 + pT,2)** — momentum sharing (subleading over sum);
  0 ≤ z ≤ 0.5, low-z distorted by the 10 GeV track-jet pT threshold.
- **Δθ_ppg,gbb** — orientation of the splitting plane relative to the gluon
  production plane (sensitive to gluon polarization).
- **log(mbb/pT)** — dimensionless pair mass from the two track-jet four-vector sum.

**Background estimation = flavour-fraction fit (§6).** Per-bin of each observable,
fit the **signed transverse impact-parameter significance s_d0 = s_j·|d0|/σ(d0)**
of tracks (sign s_j set by track vs jet axis). For each track-jet, use **s^sub_d0**
= the s_d0 of the track with the *second-largest* |d0| significance (more stable,
nearly pT-independent → inclusive fit). Binned maximum-likelihood fit of
s^sub_d0(j1) and s^sub_d0(j2); the 2D pdf factorizes into a product of marginals
(linear correlation < 5%). Three templates from MC truth labeling:
**BB** (signal: both jets b), **B** (one b + one light/c, incl. a fully-merged
g→bb̄ in one track-jet), **L+C** (rest). Fit range s^sub_d0 ∈ [−40, 70].
Background (non-BB) subtracted before unfolding. **Pythia pre-fit fractions
disagree significantly with the fitted (data) fractions**: Pythia slightly
**overestimates the BB (signal) fraction** in all cases, and the B vs L+C
fractions are inverted between Pythia and data (example bin: pre-fit BB 20% → post-
fit 17%; background 79.6% → 82.8%; χ²/dof 72/22 → 13.5/22).

**Unfolding (§7).** Iterative Bayes (RooUnfold), 4 iterations, with
acceptance/efficiency corrections for selection-level mismatches.

**Uncertainties (§8, Table 1; normalized differential cross-sections).** Dominated
by jet energy scale, unfolding, and **theoretical modeling**:

| Observable | Calo JES | Flavor tag | Tracking | Bkg fit | Unfolding | **Theory model** | Stat | **Total** |
|---|---|---|---|---|---|---|---|---|
| ΔR(b,b) | 2–3% | <1% | 1–2% | 1% | 2–3% | **3–10%** | 1% | **3–10%** |
| Δθ_ppg,gbb | 2–3% | <1% | 1–2% | 1% | 2% | **2–13%** | 1% | **3–10%** |
| z(pT) | 2–6% | <1% | 2–4% | 1–2% | 2–4% | **3–10%** | 2% | **3–14%** |
| log(mbb/pT) | 2–4% | <1% | 1–2% | 2% | 2–5% | **4–11%** | 1% | **4–12%** |

Theory-modeling uncertainty from Pythia↔Sherpa fragmentation differences (response
matrix + correction factors), added in quadrature; flavour-tagging largely
cancels because fractions are constrained in situ.

**Results (§9, Fig. 6).** **Sherpa 2.1 generally describes data better than
Pythia 8.230**, but both disagree at low log(mbb/pT), low z(pT), and across all
Δθ_ppg,gbb. The Δθ_ppg,gbb data shape is *inverted* relative to Pythia (minimum,
not maximum, at π/2); turning off Pythia gluon-polarization azimuthal asymmetry
moves it toward Sherpa/data. Pythia A14 FSR variations tested (Var2± ≈ ±10% on
final-state shower α_s; renormalization scale m²_bb/4 instead of p²_T,bb) — **no
single variation describes all the data.** Conclusion: g→bb̄ collinear
fragmentation is poorly modeled and unconstrained; particle-level spectra released
for MC tuning.

## References worth future reading   (≤3)
1. **P. Ilten, N. L. Rodd, J. Thaler, M. Williams, *Disentangling heavy flavor at
   colliders*, Phys. Rev. D 96 (2017) 054019, arXiv:1702.02947.** — SUPPORTIVE.
   New info: phenomenological method to statistically separate g→QQ̄ from other
   heavy-flavour production mechanisms at colliders. Serves: background-composition
   framing / possible template-separation methodology (`analysis_overview` §6).
2. *(b-hadron pair production at √s = 8 TeV, arXiv:1705.03374 — already queued via
   [[gluon_splitting_flavour_excitation]]; not re-listed here.)*

The dense-environment b-tagging calibration notes (ATLAS-CONF-2016-002/039 [50,51],
CMS-PAS-BTV-13-001 [52]) are CDS conference notes (anti-bot-blocked, not citeable
papers) and their close-track/merged-cluster mechanism is already covered by
[[atlas_inner_detector_tracking]] — **not worth adding.**

## Related KB docs
- [[gluon_splitting_flavour_excitation]] — background-concept HUB (FSR g→QQ̄ vs ISR flavour excitation; signal-vs-background distinction; template-stitch). **This doc is its citeable small-angle measurement anchor.**
- [[muon_source_template_fits]] — our d0 / Δp/p purity template fits; this paper's per-track s_d0 flavour-fraction fit is the jet-track analogue (method contrast, not the same observable).
- [[atlas_inner_detector_tracking]] — dense-environment / merged-cluster tracking, the mechanism behind close-by b-jet (and our close-by muon) reconstruction.
- [[run2_dimuon_backtoback_paper]] — the dimuon side of correlated HF backgrounds (charm-vs-bottom handle, both muons pT > 4 GeV).
- `Analysis/docs/analysis_overview.md` §6 — our background & template program (ground truth, external to KB).
