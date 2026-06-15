# ATLAS Run 2 Muon Reconstruction, Identification & Efficiency

**Source:** *Muon reconstruction and identification efficiency in ATLAS using the
full Run 2 pp collision data set at √s = 13 TeV*, The ATLAS Collaboration,
Eur. Phys. J. C **81** (2021) 578 — arXiv:2012.00578 / CERN-EP-2020-199.
Local PDF: `ATLAS_Run2_muon_reco_eff_EPJC81_578.pdf` (this directory).

**Scope of measurement:** 139 fb⁻¹ of pp at √s = 13 TeV (2015–2018).
Full muon acceptance |η| < 2.7, transverse momenta 3 GeV < pT < 250 GeV.
Efficiencies measured with Z→µµ (high pT) and J/ψ→µµ (low pT) tag-and-probe,
plus a double-ratio method for 2.5 < |η| < 2.7. Uncertainties reduced by ~5×
relative to the early-Run-2 publication (Ref. [7]).

> **Relevance to this analysis:** This is the *pp* muon performance paper. Our
> dimuon analysis uses pp24 (2mu4 trigger) and PbPb. The reconstruction
> strategies, working-point definitions, and tag-and-probe / scale-factor
> methodology here are the reference for understanding muon efficiency and the
> data/MC scale factors applied in pp. PbPb muon performance differs (higher
> occupancy) and is documented separately. Note "HI" in our branch names always
> means Heavy Ion — do not conflate with the pp results below.

---

## 1. Muon Reconstruction Procedure

Muon signature = minimum-ionising particle: a track in the Muon Spectrometer
(MS) and/or characteristic calorimeter energy deposits. Reconstruction relies
primarily on the Inner Detector (ID, |η| < 2.5) and the MS (|η| < 2.7), with the
calorimeters used for energy-loss corrections and for MS-independent tagging.

### 1.1 MS stand-alone track reconstruction
1. Find short straight-line **segments** in individual MS stations via a Hough
   transform (hits in MDT/CSC precision chambers + RPC/TGC trigger chambers).
2. Combine segments across stations into preliminary track candidates using a
   loose IP-pointing constraint and a **parabolic trajectory** (first-order
   approximation of magnetic bending — improved over the old straight-line
   method).
3. Combine bending-plane precision measurements with second-coordinate
   measurements from trigger chambers → 3D candidates.
4. Global χ² fit through the B-field, accounting for material interactions and
   chamber misalignments (now via **constrained nuisance parameters** for
   translational/rotational chamber displacements).
5. Outlier hit removal, re-fit, ambiguity resolution (drop tracks sharing many
   hits with higher-quality tracks; exception preserves boosted low-mass dimuon
   efficiency). Final re-fit with loose IP constraint + calorimeter energy loss,
   back-extrapolated to the beamline; pT expressed at the IP.

### 1.2 Muon types (five reconstruction strategies)
Global reconstruction combines ID + MS + calorimeter information into five types:

| Type | Abbr. | How built | Notes |
|------|-------|-----------|-------|
| **Combined** | CB | MS track matched to ID track, combined refit incl. calo energy loss | Highest quality; >98% of Medium prompt muons are CB. For \|η\|>2.5, MS tracks may combine with pixel+SCT segments → **silicon-associated forward (SiF)** muons |
| **Inside-out** | IO | ID track extrapolated outward, ≥3 loosely-aligned MS hits, combined fit | Does *not* need a standalone MS track; recovers efficiency at low pT and in limited-MS-coverage regions |
| **MS-extrapolated** | ME | Unmatched MS track extrapolated to beamline | Extends acceptance beyond ID, up to \|η\| = 2.7 |
| **Segment-tagged** | ST | ID track extrapolated to MS, tight angular match to ≥1 MS segment | Parameters from ID fit |
| **Calorimeter-tagged** | CT | ID track extrapolated through calo, energy deposit consistent with MIP | Parameters from ID fit; pT > 5 GeV required (background); retuned for \|η\| < 0.1 gap |

Most reconstruction algorithms use ID tracks down to pT = 2 GeV (CT: 5 GeV).

**Key improvements vs. Ref. [7]:** parabolic pattern recognition; SiF muons near
ID acceptance edge; alignment nuisance parameters in the fit; retuned
calo-tagging (+ a looser CT working point for tag-and-probe).

---

## 2. Identification Working Points (WPs)

A WP = a set of requirements on ID/MS hit counts, track-fit quality, and
ID–MS compatibility variables, defined per muon type. Goal: balance prompt-muon
efficiency, momentum resolution, and rejection of non-prompt muons (light-hadron
in-flight decays produce kinked, lower-quality tracks; heavy-flavour b/c decays
produce good tracks rejected mainly via vertex/isolation).

### 2.1 Key discriminating variables
- **q/p compatibility** = |q/p_ID − q/p_MS| / √(σ²(q/p_ID) + σ²(q/p_MS)) — ID/MS
  momentum agreement (CB/IO muons with MS track).
- **ρ′** = |pT,ID − pT,MS| / pT,CB — fractional ID–MS pT difference.
- **Precision station** = MS station with ≥3 MDT/CSC hits. **Precision hole
  station** = expected station with <3 hits and ≥3 missing.
- **Common ID requirement** (all WPs, CB/IO/ST/CT): ≥1 pixel hit, ≥5 SCT hits,
  ≤2 missing silicon hits. (SiF exception: ≥1 pixel, ≥4 pixel+SCT total.)

### 2.2 The three standard WPs (nested: Tight ⊂ Medium ⊂ Loose)
- **Loose** — highest efficiency, designed for H→4µ. Adds CT/ST muons in
  |η| < 0.1 gap and low-pT ST-confirmed IO muons (|η|<1.3, pT<7 GeV) on top of
  Medium. ~97% CB+IO, ~1.5% CT/ST. +20% efficiency vs Medium at 3–5 GeV.
- **Medium** — workhorse, small systematics. Accepts CB+IO with ≥2 precision
  stations (≥1 allowed in |η|<0.1). q/p compatibility < 7. Extends to
  2.5<|η|<2.7 via ME/SiF with ≥3 precision stations.
- **Tight** — highest purity, for non-prompt-limited analyses. CB+IO, ≥2
  precision stations, normalised χ² < 8, pT/|η|-dependent cuts on q/p
  compatibility and ρ′. ~6% prompt-efficiency loss but >50% background reduction
  vs Medium in 6–20 GeV.

### 2.3 Two extreme-phase-space WPs
- **High-pT** — optimal momentum resolution for pT > 100 GeV (Z′/W′ searches).
  Subset of Medium; ≥3 precision stations (≥4 in B-field inversion zones; 2-station
  allowed only if missing hits in inner station and |η|<1.3). Rejects badly-aligned
  regions (all of 1.0<|η|<1.1; support structures near φ=−1.2, −2.0). Adds a
  **σ_rel(q/p)** cut on momentum-uncertainty. Overall Z/γ* eff ≈ 80% at 100 GeV,
  ~72% at 1 TeV, ~68% at 2 TeV.
- **Low-pT** — for the softest muons (b-physics, SUSY compressed spectra). CB+IO
  only (IO must be ST-confirmed); ≥1 precision station (≥2 for |η|>1.3 with
  pT>3 GeV). Identical to Medium above pT = 18 GeV. Two variants:
  - **cut-based:** cuts on MBS (momentum-balance significance), SNS (scattering-
    neighbour significance), SCS (scattering-curvature significance), each < 3;
    plus Medium in |η|>1.55.
  - **multivariate (MVA):** gradient BDT, 8 variables (SCS, SNS, MBS, calo energy
    loss, # MS segments + directions, # missing middle-station precision hits),
    trained separately for CB/IO on tt̄ prompt vs light-hadron muons.
  - vs Medium at 3–5 GeV: +16% (cut) / +18% (MVA) prompt muons, only +0.2%/+0.1%
    light hadrons.

### 2.4 Vertex association (impact-parameter) criteria
- **d₀ significance** |d₀|/σ(d₀) < 3 (transverse, measured w.r.t. beamline; beam
  width folded into σ(d₀)).
- **|z₀ sin θ| < 0.5 mm** (longitudinal, w.r.t. primary vertex).
- Rejects in-flight hadron decays, pile-up, cosmics. Resolution asymptotes ~10 µm
  (transverse) / 50 µm (longitudinal) for pT > 10 GeV.

### 2.5 Isolation WPs
Isolation = (ΣET or ΣpT in a cone)/pT_µ. Three flavours: track-based (ID tracks),
calorimeter-based (topo-clusters), particle-flow (combination).
- **Track variables:** `ptcone20` (ΔR=0.2) or `ptvarcone30` (ΔR=min(10 GeV/pT,0.3),
  for nearby-object topologies). Min track pT 500 MeV or 1 GeV. Largely pile-up
  independent.
- **Calo:** `topoetcone20` (ΔR=0.2), pile-up corrected via ambient energy density;
  poorer resolution, more pile-up dependent.
- **Particle-flow:** track isolation + 0.4·(neutral PFlow ET in ΔR=0.2,
  `neflowisol20`); removes charged double-counting.
- **Eight isolation WPs** (Table 2): PflowLoose, PflowTight, Loose, Tight,
  HighPtTrackOnly, TightTrackOnly, PLBDTLoose, PLBDTTight. The last two use a
  **prompt-lepton BDT** (8 inputs incl. b-tagging DL1mu/RNNIP) tuned to match
  TightTrackOnly / Tight prompt efficiency with stronger HF rejection.

---

## 3. Efficiency Measurement Procedure & Results

### 3.1 Tag-and-probe (|η| < 2.5)
Select dimuon pairs (Z or J/ψ). **Tag** = stringent ID + triggered the event.
**Probe** = the test leg, reconstructed by a detector subsystem independent of the
one under study. Efficiency ε^X(p|s) = N_matched / N_all (probes matched within
ΔR < 0.05 to a muon reconstructed/identified by algorithm X). Backgrounds
subtracted via fits to the tag–probe invariant-mass spectrum.

**Probe types:** ID probes (ID tracks → test MS/ID algorithms), MS probes (ME
tracks → test ID reco), CT probes (ID + calo-tag, purer than ID probes), ST probes
(ID + segment-tag, available for pT < 5 GeV unlike CT), two-track probes (MS track
within ΔR=0.05 of an ID track → combined reco efficiency), Loose probes
(→ isolation & vertex-association efficiency).

**Scale factor:** `SF = ε_Data(X) / ε_MC(X)` — corrects simulation; method biases
cancel in the ratio. This is the quantity applied in physics analyses.

**Z→µµ (high pT):** mass window 61–121 GeV (86–96 for vertex assoc., 81–101 for
isolation). Tag: Medium, pT>27 GeV, |η|<2.5, Tight isolation, single-µ trigger.
Probe: pT>3 GeV (>10 for reco/ID). Reco/ID efficiency factorised explicitly into
four terms (Eq. 2): ε(X|ID∧MS) via two-track probes, ε(X∧¬MS|CT) and ε(MS|CT) via
CT probes, ε(ID|MS) via MS probes — residual bias ~5× smaller than Ref. [7].
Non-prompt background modelled by Eq. (3) (template with Λ ≈ 2.5× upper mass
bound, shape from same-charge pairs).

**J/ψ→µµ (low pT, 3–20 GeV):** mass window 2.7–3.5 GeV. Tag: Tight, pT>6 GeV.
ST probes (work below 5 GeV where CT fail). Simultaneous max-likelihood fit
(Crystal Ball signal + polynomial background).

**Double-ratio method (2.5 < |η| < 2.7):** no ID coverage, so no two-detector
tag-and-probe. SF from Eq. (4): data/MC ratio of Z→µµ with one forward muon
(2.5<|η|<2.7) over the same ratio with the forward muon in 2.2<|η|<2.5. Central
muon: Medium, pT>25 GeV, Tight iso. Only SFs (not absolute data efficiency) are
obtained here.

### 3.2 Systematic uncertainties (reco/ID, Z→µµ)
T&P method bias (dominant below 100 GeV; halved gen-level data/MC difference,
< 0.1% with two-track probes), probe matching, template shape (a₁,a₂ in Eq. 3),
Λ-SC (±20%), background fit, cross-section, luminosity (1.7%). Above 150 GeV the
non-prompt background estimate dominates. For J/ψ a **fit-model** systematic
dominates the lowest-pT bins. Isolation adds: Probe-PID, Mass-window, Jet-modelling
(Sherpa vs Powheg+Pythia — largest for PflowLoose), ΔR(jet,µ).

### 3.3 Key results
**Reconstruction & identification (Table 1, tt̄ MC, |η|<2.5):** prompt
efficiency ε_µ / light-hadron mis-ID ε_had —

| WP | 3–5 GeV | 5–20 GeV | 20–100 GeV | >100 GeV |
|----|---------|----------|------------|----------|
| Loose | 90 / 1.17 | 98 / 1.06 | 99 / 0.25 | 98 / 0.12 |
| Medium | 70 / 0.63 | 97 / 0.85 | 97 / 0.17 | 97 / 0.07 |
| Tight | 36 / 0.15 | 90 / 0.38 | 93 / 0.12 | 93 / 0.04 |
| Low-pT (cut) | 86 / 0.82 | 95 / 0.71 | 97 / 0.17 | 97 / 0.07 |
| Low-pT (MVA) | 88 / 0.73 | 96 / 0.66 | 97 / 0.17 | 97 / 0.07 |
| High-pT | 45 / 0.34 | 79 / 0.60 | 80 / 0.13 | 80 / 0.05 |
(values in %; Tight selects nothing below 4 GeV.)

- **Loose/Medium efficiency > 98%** for 0.1<|η|<2.5; data/MC agreement ~0.5% on
  average. **Tight > 95%**, data/MC ≤ 1%.
- For Medium SFs across (η,φ): ~0.5% agreement; only 2 bins (in |η|<0.1) deviate
  >5%. Local inefficiencies: detector support structures (|η|<1.0, φ≈−1.2, −2.0);
  misaligned MDT in 1.0<|η|<1.3; faulty CSCs in 2017/2018 (visible drops at
  specific (η,φ)); High-pT loses |η|>2.0 tracks due to strict CSC hit requirements.
- Z↔J/ψ measurements agree in the 10–20 GeV overlap. Reco/ID SFs computed for
  pT>15 GeV from Z (16 η bins following MS chamber layout); 3–15 GeV from J/ψ
  in (pT,η). Low-pT MVA gives larger forward efficiency and smaller SF
  uncertainties than cut-based; data/MC differences >10% only for |η|>2.0,
  pT<4 GeV (faulty CSCs + lower CSC segment reco in data).
- **Vertex association:** efficiency > 97% (→99% above 20 GeV); excellent
  data/MC, largest deviation ~2% at low pT near TRT edge (|η|≈1.9).
- **Isolation (Table 3):** SFs ≈ unity, per-mille uncertainties for pT>20 GeV
  away from jets; up to ~5% at low pT / near jets (Jet-modelling dominated). HF
  suppression factor (1/ε) from 8(12) at low pT to 20(100) at >25 GeV for
  PflowLoose(PflowTight); max 250 with PLBDTTight near 30 GeV. SFs delivered in
  4 ΔR(jet,µ) bins.
- **Stability:** reco/ID efficiency for Medium is flat vs integrated luminosity
  and vs ⟨µ⟩ (pile-up) across Run 2. Vertex-association efficiency shows a small
  progressive decrease (worse IP resolution at higher pile-up), stabilising mid-2017.
  Isolation stable.

### 3.4 Final precision (Conclusions)
- Z→µµ: reco/ID better than per-mille for pT>10 GeV in most regions.
- J/ψ→µµ: extends to 3 GeV, better than 1% in 5–20 GeV.
- Vertex association: better than 0.2% over full pT, < 0.01% above 20 GeV.
- Isolation (8 WPs): SFs ≈ 1, per-mille for pT>20 GeV away from jets.

---

## Detector context (Section 2)
- **ID** (|η|<2.5, 2 T solenoid): pixel (incl. IBL) + SCT + TRT (TRT to |η|=2.0).
- **MS** (air-core toroids, 2–6 T·m): three stations of **MDT** (innermost replaced
  by **CSC** for |η|>2.0). MDT/CSC single-hit resolution ~80/60 µm (bending plane).
  Trigger chambers: **RPC** (barrel, |η|<1.05) + **TGC** (endcap, 1.0<|η|<2.4),
  spatial res. 5–10 mm, also give the non-bending "second coordinate".
- Trigger: L1 hardware (<100 kHz) + HLT software (~1 kHz to disk).
- Data: 139 fb⁻¹, ⟨µ⟩=34 (13/25/38/36 for 2015/16/17/18). Z→µµ via single-muon
  triggers (pT thresh. 20→26 GeV, to avoid biasing the probe); J/ψ→µµ via
  muon+track triggers (mass 2.5–4.3 GeV, track from MS or ID, pT 3–6 GeV).
- MC: Z→µµ Powheg+Pythia8 (CT10, AZNLO), ~210 M; J/ψ→µµ Pythia8 LO (+Photos++),
  ~420 M; full Geant4 + pile-up overlay.

---

## Related KB docs   (knowledge graph)

> **This is the canonical muon-performance hub** (5 muon types, working points,
> tag-and-probe). The trigger and Run-3 docs defer their shared definitions here.

- [[atlas_run2_muon_trigger]] — companion: muon **trigger** efficiency (this doc is reco/ID); shared method = tag-and-probe.
- [[atlas_run3_muon_performance]] — Run-3 reco+trigger; reuses these muon-type/WP definitions (NSW added).
- [[atlas_inner_detector_tracking]] — ID dense-environment tracking (merged/shared clusters, ambiguity solver) underlying close-track muon reco.
- [[../../analysis/run2_hf_muon_raa]], [[../../analysis/run2_hf_muon_vn]], [[../../analysis/run2_dimuon_backtoback_paper]] — analyses that apply this tag-and-probe efficiency method.
- [[../../concepts/muon_source_template_fits]] — uses the d0/z0 resolution and ID/MS momentum measurements described here.
