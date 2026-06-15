# ATLAS Inner-Detector Tracking — Dense Environments & High-Occupancy

**Primary source (this doc):** *Performance of the ATLAS track reconstruction
algorithms in dense environments in LHC Run 2*, The ATLAS Collaboration,
Eur. Phys. J. C **77** (2017) 673 — CERN-EP-2017-045, DOI
10.1140/epjc/s10052-017-5225-7.
**URL:** https://arxiv.org/abs/1704.07983 (PDF https://arxiv.org/pdf/1704.07983)
**Classification:** SUPPORTIVE (no local PDF committed — kept lean; URL above).
**Added:** 2026-06-15

> This doc covers ID **tracking** physics that the pp muon-reco PRIMARY doc
> ([[ATLAS_Run2_muon_reconstruction]]) does NOT explain: (1) the hit-level
> merged/shared-cluster + ambiguity-solving mechanism by which two **close-by**
> tracks have **correlated** reconstruction (relevant to our pair ε_reco and its
> ΔR axis), and (2) a pointer to ID tracking in the **high-occupancy Pb+Pb**
> environment (relevant to the centrality dependence of our PbPb reco
> efficiency). Both are **[background-for-writing]** — we do **not** import any
> numbers from these sources; our reco efficiency is measured in our own
> fullsim/overlay MC.

---

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Merged vs shared clusters; ambiguity solver limits shared clusters (≤2 tracks/cluster, ≤1 shared pixel cluster) → close tracks compete for hits | Physical justification that ε_reco is a single **pair** efficiency binned in (pair pT, pair η, **ΔR**), **not** ε₁·ε₂ — `analysis_overview.md` §4b ("shared ID hits, ambiguity resolution… the ΔR axis captures that correlation"); roadmap step 12; IntNote reco-eff section / thesis | [background-for-writing] |
| ID layout & clusterization (pixel+IBL+SCT+TRT, CCA, ToT/dE/dx) at the hit level | Context for HIJING-**overlay** reco (overlay injects HIJING at the **hit** level) and for why merged clusters arise | [background-for-writing] |
| (separate slice, below) ID charged-particle tracking efficiency in **high-occupancy Pb+Pb** is lower and **centrality-dependent** | Motivates the centrality dependence of our PbPb pair reco efficiency — `analysis_overview.md` §4b ("central Pb+Pb… genuinely lower… real occupancy effect from the underlying event"); thesis/IntNote | [background-for-writing] |

The single most analysis-specific point: this source documents *why* close-by
tracks lose efficiency at the hit level, which is the documented reason our
reco-efficiency is a **pair** efficiency with a ΔR axis (see the project memory
"pair_reco_effcy_definition" and `analysis_overview.md` §4b).

## Scope & condition-difference warnings

Primary source: **Run 2 pp at √s = 13 TeV**, 3.2 fb⁻¹, IBL installed; the
dense-environment regime studied is the **cores of 200–1600 GeV jets**, where
charged-particle separations approach the **sensor granularity** (a few pixels).
- Our signal is a low-mass b→μμ pair with **pair pT > 8 GeV** and **ΔR > 0.05** —
  a far **less dense / lower-pT** regime than TeV-scale jet cores. The *mechanism*
  (shared/merged clusters, ambiguity solving) is the same physics; the
  *magnitude* of any close-track inefficiency does NOT transfer. ACKNOWLEDGE the
  mechanism; do NOT assume the size/direction of any efficiency loss for our
  pairs (KB GUIDE §5).
- Run 2 detector/reco config; our analysis is **Run 3 5.36 TeV** with HIJING
  overlay — conditions differ; do not quantify the difference.
- HI charged-particle tracking efficiencies (slice below) are for **primary
  charged hadrons**, NOT for our **muon-pair** reconstruction — different object,
  different selection. We measure our own pair efficiency in overlay MC.

## Content summary

### 1. ID at the hit level (Sec. 2–3.1)
ID covers |η| < 2.5: silicon **pixel** (incl. **IBL**, mean r ≈ 33 mm, pixel
50 µm × 250 µm; B-layer at r ≈ 50.5 mm; further layers 88.5, 122.5 mm) + **SCT**
(4 double strip layers, r 299–514 mm, 80 µm pitch, ~4 3D measurements/track) +
**TRT** (to r ≈ 1082 mm, |η| < 2.0). Pixels read out **charge via Time-over-
Threshold (ToT)**, proportional to deposited energy → enables a **dE/dx** measure.
**Clusterization:** a connected-component analysis groups above-threshold
pixels/strips into clusters → 3D space-points (one per pixel cluster; SCT combines
both sides). In **dense environments** the spatial separation between two charged
particles is only a few pixels, so CCA sometimes makes **one cluster holding
deposits from multiple particles**.

### 2. Cluster classes (Sec. 3.1) — the key vocabulary
- **single-particle cluster** — charge from one particle (truth).
- **merged cluster** — charge from multiple particles (truth).
- **identified as merged** — reco flags a cluster as compatible with merging
  (via a neural-network cluster splitter).
- **shared cluster** — used in ≥2 reconstructed tracks but NOT identified as
  merged. Multiply-used clusters are *either* identified-as-merged *or* shared.

### 3. Iterative track finding + ambiguity solving (Sec. 3.2–3.3) — the mechanism
Seeds from 3 space-points → combinatorial **Kalman filter** builds candidates.
Candidates are then rated by a **track score** (clusters add score per
subdetector resolution; **holes** subtract; large χ² penalized; log(pT) promotes
energetic tracks) and processed in **descending score** in an **ambiguity
solver**. Handling of multiply-used clusters — this is the close-track
correlation mechanism:
- A candidate is compared only to **already-accepted** tracks.
- A cluster may be **shared by at most two tracks**, preference to the
  higher-scored (earlier-processed) track.
- A track may have **at most two shared clusters**.
- Quality cuts to be accepted: pT > 400 MeV, |η| < 2.5, ≥7 pixel+SCT clusters
  (12 expected), **≤1 shared pixel cluster OR ≤2 shared SCT clusters on the same
  layer**, ≤2 holes.
Consequence: when two muons are close enough to **share** silicon clusters, the
ambiguity solver can strip clusters from / reject the lower-scored track → the
two-muon reconstruction is **correlated**, not independent. (For isolated
primaries reco eff is very high; the paper notes muon reco eff > 99% for isolated
tracks, citing the Run 1 muon paper.)

### 4. Data-driven inefficiency in jet cores (Sec. 6) — magnitude in the DENSE regime only
A novel data-driven method uses pixel **dE/dx** to tag clusters from two charged
particles. The measured **fraction of charged particles (creating such clusters)
that fail to be reconstructed**: **0.061 ± 0.006 (stat) ± 0.014 (syst)** for jet
pT 200–400 GeV and **0.093 ± 0.017 (stat) ± 0.021 (syst)** for jet pT
1400–1600 GeV (Abstract; rises with jet pT). **These numbers are for TeV-scale
jet cores and do NOT apply to our b→μμ pairs** — quoted only to show the
mechanism is real and pT/density-dependent.

### 5. High-occupancy Pb+Pb tracking — RELATED SLICE (different sources; background only)
The dense-environment paper is **pp jet-core** dense, not Pb+Pb **occupancy**
dense; these are distinct regimes (collimated few-track vs. global high
multiplicity). ATLAS HI charged-particle analyses determine ID **track-
reconstruction efficiency vs (pT, η, centrality)** and a **fake-track rate**, and
report that efficiency is **lower in central than peripheral** collisions and at
**forward vs mid-rapidity**, with fakes largest at low pT / high centrality. This
is the documented ID-level basis for the **centrality dependence** of tracking in
Pb+Pb. **Specific per-centrality efficiency numbers are NOT verified here** —
see Future-Read; treat any number as ⟨unverified — check the source PDF⟩.
Representative sources (HI charged-particle spectra / flow performance sections):
- arXiv:1808.03951 (charged-particle anisotropy, 5.02 TeV Pb+Pb).
- arXiv:2412.15658 (charged-particle anisotropy at high pT, 5.02 TeV Pb+Pb).
Again: these are **charged-hadron** tracking efficiencies, NOT our muon-pair reco
efficiency, which we measure in HIJING overlay MC.

## References worth future reading   (≤3)
1. **ATLAS HI charged-particle tracking performance** (a dedicated Pb+Pb
   tracking-efficiency-vs-centrality source, e.g. the performance section of
   arXiv:1808.03951 or 2412.15658) — SUPPORTIVE. *New info:* verified
   centrality-dependent ID track efficiency & fake rate numbers in the
   high-occupancy regime. *Serves:* writing the PbPb reco-eff motivation
   (thesis/IntNote); cross-check that our overlay efficiency trend with
   centrality is qualitatively expected.
2. **ATLAS HF semileptonic muons in Pb+Pb** (arXiv:1805.05220, 2.76 TeV;
   arXiv:1209.6282) — PRIMARY (different topic: prior measurement, not ID
   tracking). *New info:* HI-optimized muon hit selections and a prior
   single-b/HF muon suppression result to **compare against**. *Serves:* Intro /
   results comparison — belongs in an *analysis-papers* KB doc, not this tracking
   doc.
3. None further worth adding from this source's reference list (the rest is
   generic ID/jet-reco machinery already covered by [[ATLAS_Run2_muon_reconstruction]]).

## Related KB docs
- [[ATLAS_Run2_muon_reconstruction]] — pp muon reco/ID; covers ID **track-quality
  cuts**, **d0/z0** significance + resolution, vertex association, and one bare
  sentence on the ambiguity-resolution exception for boosted low-mass dimuons.
  **This doc supplies the missing hit-level mechanism** behind that sentence and
  behind our pair (ΔR-binned) ε_reco.
