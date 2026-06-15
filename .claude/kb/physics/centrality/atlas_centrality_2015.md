# ATLAS Centrality Determination — 2015 Pb+Pb (5.02 TeV)

**Source:** D. Perepelitsa, P. Steinberg, M. Zhou, S. Radhakrishnan, S. Tapia Araya,
"Centrality determination in √s_NN = 5.02 TeV Pb+Pb collision data in 2015" —
ATLAS internal supporting note, ATL-COM-PHYS-2016-XXX, Draft v0.1, 4 Sept 2016.
**arXiv / DOI:** none (internal draft note).
**PDF:** `./centrality 2015.pdf` (filename has a space — quote it)
**Classification:** PRIMARY (foundational ATLAS FCal-centrality methodology; 18 pp)
**Added:** 2026-06-15

This is the **earliest detailed ATLAS centrality supporting note** in the KB and the
**canonical methodology reference** for how ATLAS turns FCal ΣE_T into centrality and
Glauber geometry (⟨N_part⟩, ⟨N_coll⟩, ⟨T_AA⟩). It is the **methodology-depth companion**
to the centrality *values* paper — see [[atlas_centrality_2023]] for the numbers our
analysis actually applies. **All numbers here are 5.02 TeV / 2015 (Run 2);** the method
transfers, the numbers do not (see Scope warning).

## Relevance to this analysis   (specific)
| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| FCal ΣE_T as the centrality estimator (3.1<\|η\|<4.9); event/pileup selection for the MB ΣE_T spectrum | Centrality determination — overview §4e, roadmap step 6, IntNote §6 | [method-we-use] |
| Glauber two-component model (TCM) fit to the ΣE_T distribution; working-point + efficiency; centile (ΣE_T-cut) extraction | The procedure behind our 2023-Glauber centrality thresholds & cross-year FCal scaling | [method-we-use] |
| Extraction of ⟨N_part⟩, ⟨N_coll⟩, ⟨T_AA⟩ = ⟨N_coll⟩/σ_NN per centile; their systematic-uncertainty recipe | ⟨T_AA⟩ enters R_AA (overview §3, §4d); how its uncertainty is built | [method-we-use] |
| Why centrality is defined only above a ΣE_T working point (removes photonuclear/diffractive non-Glauber events) | Justifies the analysis MB/event-selection framing in the note | [background-for-writing] |

## Scope & condition-difference warnings
**System:** Pb+Pb, **√s_NN = 5.02 TeV**, **2015** data (Run 2), 35 MB runs (33 after GRL),
~0.52 nb⁻¹. Our analysis is **Run 3, 5.36 TeV**, using **2023 Glauber thresholds**
([[atlas_centrality_2023]]). ACKNOWLEDGE the era/energy difference; **DO NOT** carry any
ΣE_T cut, fit parameter, x value, efficiency, or ⟨T_AA⟩/⟨N_part⟩/⟨N_coll⟩ number from this
note into the 5.36 TeV analysis — those are tabulated for 5.02 TeV/2015 only. What transfers
is the **method**. Do not estimate the size/direction of the energy/era difference (§5).

## Content summary

### 1. Centrality estimator & event selection  [method-we-use]
Centrality is characterized by **ΣE_T summed over both FCal modules at 3.1 < |η| < 4.9**
(same observable as Run 1). The MB ΣE_T distribution (p.2–3) is built from a logical-OR of
two MB triggers (`HLT_noalg_mb_L1TE50`, `HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50`), prescale-weighted,
requiring GRL + a reconstructed vertex.
- **Pileup rejection (§2.1):** uses the **FCal-energy vs ZDC-neutron anti-correlation** —
  pileup (two overlaid Pb+Pb) forms a separate band. `HIPileupTool` cuts are tuned so
  ≤0.1% of good events are rejected at any energy; >80% pileup rejection above 0.5 TeV,
  ≈100% above 4.5 TeV (very central). (Mechanism is the same physics our PbPb event
  selection's FCal/ZDC and preamp cuts rest on — but our cut values are independent.)

### 2. Centrality (Glauber two-component) model  [method-we-use]
- **Glauber MC:** 10⁶ PHOBOS Glauber MC events at 5.02 TeV, **σ_NN = 70 mb**, Woods-Saxon
  **(R, a) = (6.62 fm, 0.546 fm)**, nucleon hard core **d_N = 0.4 fm** → per-event N_part, N_coll (p.3–4).
- **Two-component model (TCM):** soft particle production (here ΣE_T at large η) scales with
  the **number of "fundamental elements"** (eq. 1):

  **N_elem = x·N_coll + (1−x)·N_part/2**

  x weights the binary-collision vs participant scaling; for pp (N_part=2, N_coll=1) N_elem=1
  for any x, so one element ≈ one pp collision's ΣE_T at the same √s.
- **ΣE_T response per element** (eqs. 2–5): the ΣE_T from a single element is modeled as a
  **Gaussian (default)** P(ΣE_T) = Gaus(ΣE_T; σ_pp, µ_pp), or a **Gamma (alternate)**
  Gamma(ΣE_T; θ_pp, k_pp). An event with N_elem elements is the **N_elem-fold convolution**
  (Gaussian scales µ→µ·N_elem, σ→σ·√N_elem). The full PbPb ΣE_T distribution (eq. 6) is the
  sum over N_elem of the N_elem-fold response weighted by P(N_elem) from Glauber. This model
  shape is **fit to the data**, floating the pp-like (µ,σ or θ,k) parameters and x.

### 3. Analysis: working point, x, fit, centile cuts  [method-we-use]
- **Working point (§4.1):** Glauber describes only the **hadronic** cross-section; data also
  contain photonuclear/diffractive events at small ΣE_T. A **reconstructed rapidity-gap
  analysis** (49 η-rings of δη=0.2 over |η|<4.9; edge-gap Δη_edge = largest run of empty rings
  from either edge; event = diffractive candidate if **Δη_edge > 2.0**) shows that above
  **ΣE_T = 40 GeV** the large-gap fraction is negligible (<1%). → **WP = 40 GeV**; the fit and
  all centrality definitions are done **only above 40 GeV**.
- **Determining x (§4.2):** two external constraints. (a) pp-scaling of ⟨ΣE_T⟩ (eqs. 7–10,
  using ⟨N_elem⟩) gives x ≈ 0.11; (b) flattening the sim/data ratio in the shoulder gives
  x ≈ 0.07. Nominal **x = 0.09 ± 0.02** (Run 1 found 0.088 by other methods).
- **Nominal fit:** fix x=0.09, float Gaussian µ,σ, fit above 40 GeV, normalize sim & data to
  equal integral above WP. Result (Fig. 5): **µ = 11.8 GeV, σ = 5.65 GeV**, and the
  **fraction of the Glauber distribution above the WP = 84.9% → rounded to 85% ("efficiency")**.
  Gamma alternate: k=1.74, θ=6.69 GeV, f=84.6%. Refitting at x=0.07/0.11 shifts f by ±1%.
- **ΣE_T centile cuts (§4.3):** MB data above the WP are sorted by ΣE_T and divided into
  **85 equal-population bins** (= the 85% efficiency); the bin edges are the percentile cuts.
- **Geometric quantities (§4.4):** after the fit, the 2-D (ΣE_T, N_part) and (ΣE_T, N_coll)
  distributions are populated with the best-fit model; the **mean N_part, N_coll per ΣE_T
  (centrality) range** are read off. **⟨T_AA⟩ = ⟨N_coll⟩ / σ_NN** (σ_NN = 70 mb except in
  two Glauber variations).

### 4. Systematic uncertainties (§5)  [method-we-use]
Two dominant sources; for each parameter (N_part, N_coll, T_AA) each up/down pair is
symmetrized (half-difference) and the sources added in quadrature.
- **Geometric (§5.1):** six alternate Glauber sets — σ_NN = 65 / 75 mb; WS (R,a) =
  (6.68, 0.536) / (6.56, 0.556) fm; hard core d_N = 0.2 / 0.6 fm. The fit is repeated with each
  alternate N_elem distribution. (Fig. 7: per-quantity ratios vs nominal in 10% bins, all ≲ a few %.)
- **Efficiency (§5.2):** assume **84% or 86%** (instead of 85%) of the Glauber distribution lies
  above the WP → a nominal centrality range maps to one scaled by 84/85 or 86/85. For analyses
  not using geometry (e.g. correlations) an alternate ΣE_T-cut set divides into 84/86 bins.
- **Negligible (§5.3):** Gaussian-vs-Gamma response (Gaussian fits data better; geometry
  differs negligibly); WP varied 30↔60 GeV — results consistent.

### 5. Results (§6) — 5.02 TeV / 2015 ONLY (do not transfer)
- **Table 1:** ΣE_T cuts for 1% centile boundaries, nominal (85%) + 84%/86% variants. E.g.
  80% cut = 0.0637 TeV, 50% = 0.525 TeV, 10% = 2.989 TeV, 1% = 4.263 TeV.
- **Table 2 (N_part), Table 3 (N_coll), Table 4 (T_AA), Table 5 (R_coll = N_coll ratios).**
  Example 10%-bin values (5.02 TeV): ⟨N_part⟩ = 358.8 / 264.1 / 189.2 / 131.4 / 87.0 / 53.9 /
  30.6 / 15.4 for 0-10…70-80%; ⟨T_AA⟩ [mb⁻¹] = 23.35 / 14.33 / 8.638 / 4.946 / 2.634 / 1.281 /
  0.565 / 0.223 for the same bins (δT_AA/T_AA ≈ 0.9% central → ~9% at 70-80%). **These are the
  5.02 TeV numbers; our analysis uses the 5.36 TeV/2023 values in [[atlas_centrality_2023]].**
- **Flatness/run-stability (§6.2):** centrality-bin event fractions checked flat vs run;
  a few caveated runs noted (early run missing VTE50 in 75–80%; two GRL-excluded runs).

## References worth future reading   (§6; ≤3)
- Loizides, Nagle, Steinberg, "Improved version of the PHOBOS Glauber Monte Carlo,"
  SoftwareX 1-2 (2015) 13, **arXiv:1408.2549** — SUPPORTIVE. The actual Glauber MC code that
  generates N_part/N_coll/T_AA used above; relevant to our T_AA provenance. **Likely already
  covered by [[glauber_modeling]]** — add only if that doc does not treat the code itself.
- (Run-1 predecessor ATL-COM-PHYS-2011-427 and the p+Pb gap note ATL-COM-PHYS-2013-588 are
  method ancestors only — lower priority; not added.)

No further references worth adding.

## Related KB docs
- [[atlas_centrality_2023]] — the centrality **values** paper; provides the 5.36 TeV / 2023
  Glauber numbers our analysis applies. This doc supplies the method depth behind it.
- [[glauber_modeling]] — Glauber MC / Woods-Saxon / N_part-N_coll-T_AA concept hub; this note
  is the ATLAS application of that model to FCal ΣE_T.
- [[hi_big_picture]] — Glauber centrality in the broader QGP/HIC picture (Intro/thesis framing).
</content>
</invoke>
