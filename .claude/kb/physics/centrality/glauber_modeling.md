# Glauber Modeling in High-Energy Nuclear Collisions

**Source:** M.L. Miller, K. Reygers, S.J. Sanders, P. Steinberg, "Glauber Modeling in High Energy Nuclear Collisions" — Ann. Rev. Nucl. Part. Sci. **57** (2007) 205.
**arXiv / DOI:** arXiv:nucl-ex/0701025 — doi:10.1146/annurev.nucl.57.090506.123020
**PDF:** `./nucl-ex_0701025.pdf` (PRIMARY)
**Classification:** PRIMARY — the foundational, canonical methodology review for the Glauber model. **This is the designated hub for Glauber methodology**: other docs (the ATLAS centrality papers, the big-picture review) defer the *mechanism* here and keep only their own LHC-specific numbers.
**Added:** 2026-06-15

The authoritative pedagogical reference for how heavy-ion experiments turn a
measured event-activity distribution into **centrality classes** and the
**Glauber geometric quantities** ⟨N_part⟩, ⟨N_coll⟩, ⟨T_AB⟩. It is the method our
centrality determination and R_AA normalization rest on. Its **method** is what
we use (`[method-we-use]`); its **numbers** are RHIC-era (Au+Au/Cu+Cu/d+Au,
√s_NN ≤ 200 GeV) and must **not** be ported to our 5.36 TeV Pb+Pb (see warnings).

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Definitions of **N_part** (participants/wounded nucleons), **N_coll** (binary NN collisions), **T_AB** (nuclear overlap function) and their relations | Centrality & R_AA normalization vocabulary (analysis_overview §4d–e); IntNote §6 centrality section | [method-we-use] |
| **Centrality from percentiles of a measured quantity** (impact parameter b is NOT measured): bin the measured event-activity distribution by fraction of the total, map the same percentiles onto a Glauber-MC distribution to read off ⟨N_part⟩/⟨N_coll⟩ per class | The principle behind our FCal ΣE_T → 2023 Glauber-threshold centrality (analysis_overview §4e) | [method-we-use] |
| **⟨T_AB⟩ = ⟨N_coll⟩ / σ_NN^inel** (Eq. 16/26) and **N_coll = AB·T̂_AB·σ_NN^inel** (Eq. 7); T_AB-scaling of hard processes (Eq. 20–24) | Definition of ⟨T_AA⟩ used in our R_AA normalization (analysis_overview §3, roadmap Q2.3) | [method-we-use] |
| **R_AB definition** R_AB = (1/⟨T_AB⟩)·(dN^AB/dp_T)/(dσ^pp/dp_T) (Eq. 27); T_AB-scaling test via direct photons (R≈1) | The normalization form of our primary observable R_AA (analysis_overview §3) | [method-we-use] |
| **Inputs:** Woods-Saxon nuclear density (Eq. 1) + measured inelastic σ_NN | What a Glauber calc needs; what the systematics vary | [method-we-use] |
| **Systematic-uncertainty method** for N_part/N_coll/T_AB (vary σ_NN, Woods-Saxon params, nucleon hard core, NN-overlap shape, trigger/centrality-cut efficiency) and the σ_NN cancellation in T_AB | How to assign the ⟨T_AA⟩ / centrality systematic in OUR note (roadmap Q2.3 needs uncertainties) | [method-we-use] |
| Optical-limit vs Glauber-Monte-Carlo: agree for ⟨N_part⟩/⟨N_coll⟩ vs b, differ in peripheral classes and fluctuations | Background for choosing/citing a Glauber implementation; IntNote/thesis | [background-for-writing] |
| History of the Glauber model; geometric scaling N_coll ∝ N_part^{4/3} | Thesis ch. background on centrality | [background-for-writing] |

## Scope & condition-difference warnings
RHIC era: Au+Au, Cu+Cu, d+Au at **√s_NN = 19.6–200 GeV**. All numeric inputs and
results here are RHIC-specific:
- **σ_NN^inel ≈ 32–42 mb** over √s_NN = 20–200 GeV (Fig. 2). Our Pb+Pb is
  **5.36 TeV** where σ_NN^inel is **larger** — do NOT use 42 mb. The source gives
  only the *method*; the 5.36 TeV value comes from elsewhere.
- Woods-Saxon parameters quoted are for **¹⁹⁷Au** (R=6.38 fm, a=0.535 fm) and
  **⁶³Cu** (Eq. 1, p.6) — NOT ²⁰⁸Pb. Use Pb parameters for our system.
- Systematic-uncertainty *sizes* (N_part ~20% peripheral → ~3% central; N_coll
  ~10% for central; σ_geo ~10%) are RHIC numbers; ACKNOWLEDGE them as the method's
  typical scale, but do NOT assume the same magnitude/direction for our 5.36 TeV
  Pb+Pb with FCal-based centrality (§5 no-assumption rule).
- Centrality here is from **mid/forward charged multiplicity or ZDC spectators**;
  OUR analysis uses **forward calorimeter ΣE_T (FCal)** — same *principle*
  (percentiles of a monotonic activity variable), different detector
  (analysis_overview §4e). The LHC/ATLAS-specific implementation lives in the
  sibling ATLAS centrality docs.

## Content summary

### Geometric quantities (§1–2; the definitions)
A nucleus–nucleus (A+B) collision is modeled as an **independent sequence of
binary nucleon–nucleon collisions** along straight-line trajectories.
- **Participant / wounded nucleon:** a nucleon that undergoes **≥1 inelastic NN
  collision**. **N_part** = number of such nucleons in A and B.
- **Spectators:** nucleons that miss (measured at beam rapidity, e.g. ZDC neutrons).
- **N_coll:** total number of binary NN collisions in the event.
- Impact parameter **b is the fundamental geometric variable but is NOT
  measurable** (femtoscopic scale); it is inferred statistically.

Inputs (§2.2): (1) **nuclear density** ρ(r), parameterized as a 3-parameter Fermi
(Woods-Saxon) distribution
```
ρ(r) = ρ0 · (1 + w (r/R)^2) / (1 + exp[(r − R)/a])
```
(Eq. 1; ρ0 central density, R radius, a skin depth, w deformation). (2) the
**measured inelastic NN cross section σ_NN^inel** — the only nontrivial
beam-energy dependence in a Glauber calc (low-q², not pQCD-calculable).

### Optical limit vs Glauber Monte Carlo (§2.3–2.5)
- **Optical limit:** treats each nucleus as a smooth continuous density; the
  thickness function for nucleus A is T̂_A(s)=∫ρ̂_A(s,z)dz (probability per unit
  transverse area). The **nuclear overlap (thickness) function** is
  ```
  T̂_AB(b) = ∫ T̂_A(s) T̂_B(s − b) d²s          (Eq. 3)   [units: 1/area]
  ```
  interpreted as the effective overlap area for an A–B nucleon pair to interact;
  interaction probability = T̂_AB(b)·σ_NN^inel. From the binomial collision
  distribution (Eq. 4):
  ```
  N_coll(b) = A·B · T̂_AB(b) · σ_NN^inel        (Eq. 7)
  ```
  and N_part(b) from Eq. 8 (each nucleon weighted by its probability of ≥1 hit).
  Total inelastic cross section from Eq. 6.
- **Glauber Monte Carlo (GMC):** place A and B nucleons stochastically from ρ(r),
  draw b from dσ/db = 2πb, count an NN collision when transverse separation
  d ≤ √(σ_NN^inel/π) (black-disk; Gaussian/gray-disk alternatives exist).
  Directly yields event-by-event N_part, N_coll, and can be coupled to a
  multiplicity model + detector simulation for an apples-to-apples match to data.
- The two **agree** for ⟨N_part⟩ and ⟨N_coll⟩ vs b out to b ≈ 2R (Fig. 5), but
  the optical limit misses **event-by-event density fluctuations** and "shadowing"
  → larger total σ, and differs from GMC in **peripheral** centrality classes
  (Fig. 16) and for fluctuation observables. GMC is the modern default.
- Geometric scaling: **N_coll ∝ N_part^{4/3}** independent of system size (Eq. 13).

### Relating Glauber to data — centrality (§3, the core method)
**Neither N_part nor N_coll is directly measured.** The mapping (§3.1):
1. Measure a per-event activity distribution dN_evt/dN_ch (charged multiplicity,
   or forward energy) that is **monotonically related to b**: central → high
   multiplicity / few spectators; peripheral → low multiplicity / many spectators.
2. **Define centrality classes as fractions of the total integral** of that
   distribution (integrate from high activity downward). E.g. the 10–20% class is
   bounded by n₁₀, n₂₀ with ∫_{n}^{∞}/∫_{0}^{∞} = 0.1 and 0.2 (Eq. 14).
3. Apply the **same percentile binning** to a simulated distribution (GMC ×
   multiplicity model, with detector response), and read off the **mean Glauber
   quantity** ⟨N_part⟩, ⟨N_coll⟩ for the simulated events in each class.
   - Robust to an overall scale mismatch between measured and simulated activity
     (the percentile boundaries need not coincide in raw units).
4. **Leading systematic = uncertainty on the total measured cross section** (i.e.
   what fraction of σ_inel a given cut actually selects): propagated by varying
   the denominator of Eq. 14 and re-extracting. Grows toward peripheral classes.
   (Example quoted: STAR Au+Au 200 GeV 10–20% ⟨N_part⟩≈234±6 from a 5% σ_tot
   uncertainty alone.)

Caveats: **auto-correlation/acceptance bias** if the same multiplicity is used
both to define centrality and as the observable per participant (§3.3) — mitigated
by using a centrality detector in a different η region from the measurement.

### Geometric systematics for N_part / N_coll / T_AB (§3.4.2)
Total systematic = quadrature sum of variations of: σ_NN^inel (±3 mb),
Woods-Saxon R and a, **nucleon hard core** (min. 0.8 fm separation), **NN overlap
function shape** (black vs gray disk vs Gaussian), centrality-cut origin, and
**minimum-bias trigger efficiency** (which percentile is actually selected). RHIC
result (Fig. 19): σ(N_part) ~20% peripheral → ~3% central; σ(N_coll) ~10% for
central. **Key for T_AB:** because
```
⟨T_AB⟩_f = ⟨N_coll⟩_f / σ_NN^inel      (Eq. 16 / 26)
```
the **σ_NN^inel uncertainty cancels** in ⟨T_AB⟩ — so ⟨T_AB⟩ carries the **same
systematics as ⟨N_coll⟩ except** the σ_NN piece.

### Hard-process (T_AB) scaling and R_AB (§4.2)
Hard processes (cross section σ_hard^pp) scale with the overlap function via the
pQCD factorization theorem: yield per A+B encounter = T_AB(b)·σ_hard^pp (Eq. 20),
T_AB normalized so ∫T_AB d²b = AB. For a centrality class f the spectrum-level
relation is (1/N_inel^AB) dN_x^{A+B}/dp_T = ⟨T_AB⟩_f · dσ_x^{pp}/dp_T (Eq. 24).
The **nuclear modification factor**:
```
R_AB(p_T) = [ (1/N_inel^AB) dN_x^{A+B}/dp_T ] / [ ⟨T_AB⟩_f · dσ_x^{pp}/dp_T ]   (Eq. 27)
```
R_AB = 1 in the absence of nuclear effects. Validated at RHIC: **direct photons
follow T_AB scaling (R≈1)** while neutral pions are suppressed in central
collisions — the experimental proof of T_AB scaling of hard probes and of
jet quenching (Fig. 23). This is exactly the ⟨T_AA⟩-normalized R_AA form our
analysis uses (analysis_overview §3).

## References worth future reading   (§6; ≤3)
**None worth adding.** The references here are 1950s–1990s Glauber-theory
foundations (e.g. Bialas–Bleszynski–Czyż wounded-nucleon model, ref 15) and
RHIC-experiment centrality papers — superseded for OUR purposes by this review
itself (the methodology) and by the **LHC/ATLAS-specific implementations** that
are already sibling docs ([[atlas_centrality_2023]], [[atlas_centrality_2015]]).
The 5.36 TeV σ_NN^inel and Pb Woods-Saxon inputs we actually need come from those
ATLAS papers, not from this RHIC-era reference list.

## Related KB docs   (knowledge graph)
- [[atlas_centrality_2023]] — **spoke.** ATLAS Run-3 Pb+Pb centrality calibration; the LHC/FCal-specific implementation of this method (our 2023 thresholds, ⟨T_AA⟩, σ_PbPb at 5.36 TeV). Defers Glauber mechanism here.
- [[atlas_centrality_2015]] — **spoke.** ATLAS Run-2 Pb+Pb centrality; earlier LHC implementation. Defers Glauber mechanism here.
- [[hi_big_picture]] — its Glauber paragraph (§3, N_part/N_coll/spectators/percentile centrality) is the one-paragraph framing version; **defers the full methodology to this doc**.
- [[run2_hf_muon_raa]] — the ATLAS HF-muon R_AA analysis ours derives from; uses this T_AA-normalized R_AA and FCal-Glauber centrality.
- [[analysis/overview]] — our R_AA observable and centrality/normalization (§4d–e) rest on this method.
