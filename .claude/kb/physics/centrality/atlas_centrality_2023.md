# ATLAS Centrality Determination — 5.36 TeV Pb+Pb, 2023 data

**Source:** F. Pauwels, M. Spousta (Charles U., Prague), D. Perepelitsa (U. Colorado, Boulder),
*Centrality Determination in √s_NN = 5.36 TeV Pb+Pb Collisions from 2023 data running* —
ATLAS internal supporting note, **ATL-COM / GROUP-2024-033**, draft v0.3, 19 Feb 2025.
**arXiv / DOI:** none (ATLAS-internal note, not public)
**PDF:** `./centrality 2023.pdf` (PRIMARY — note the space in the filename)
**Classification:** PRIMARY
**Added:** 2026-06-15

> **This is the direct source of the centrality inputs our analysis uses.** Its
> Table 7 supplies the exact ⟨T_AA⟩ values hardcoded in the code (confirmed
> numerically below), and its Glauber methodology determines the validity of our
> Pb+Pb centrality determination. All values here are for **5.36 TeV, 2023 data**.

## Relevance to this analysis   (REQUIRED — specific)
| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| **Table 7 ⟨T_AA⟩ per centrality** | R_AA normalization (analysis_overview §3 R_AA, §4d), code `PbPbBaseClass.h::make_crossx_factors_pbpb_*`, roadmap **Q2.3** | [method-we-use] |
| **Table 1 FCal ΣE_T percentile (centile) cuts** | Centrality determination from FCal ΣE_T (analysis_overview §4e "2023 Glauber thresholds"), roadmap step 6 | [method-we-use] |
| **Glauber-MC parameters** (σ_NN, Woods-Saxon, hard core, TCM x) | Provenance/validity of the centrality calibration; IntNote §6 centrality, systematics | [method-we-use] |
| Tables 5/6 ⟨N_coll⟩, ⟨N_part⟩ + geometric/efficiency systematics | IntNote centrality table; thesis/note geometry quantities | [method-we-use] |
| Total Glauber-like σ above WP / efficiency (89%) | Centrality definition & MB normalization context | [background-for-writing] |

## Scope & condition-difference warnings
- **System / energy:** Pb+Pb, **√s_NN = 5.36 TeV**, ATLAS **2023** data running.
  *(The intro §1 once writes "5.32 TeV" — a typo in this draft; the title, the
  Glauber simulation §3, and all systematics use 5.36 TeV.)*
- **Our analysis is also 5.36 TeV (Pb+Pb 2023/2024/2025)** → energies MATCH for
  2023. **Caveat (roadmap Q2.3):** there is **no official 2024/2025 centrality
  determination yet** — our analysis **reuses these 2023 T_AA as a placeholder**
  for 2024/2025 (cross-year FCal scaling brings 24/25 ΣE_T onto the 2023 scale,
  analysis_overview §4e). State this explicitly in the note.
- **Draft status:** v0.3, "ATLAS DRAFT". Numbers may shift in a later approved
  version; treat as the current best internal reference.

## Content summary

### Method (FCal ΣE_T + Glauber, §1–4)
Centrality = sum of transverse energy in **both FCal modules at 3.2 < |η| < 4.9**
(ΣE_T). Follows the Run 1/Run 2 procedure but redone for 5.36 TeV / 2023
conditions. Data: PC + CC streams, nominal MinBias triggers
(`HLT_mb_sptrk_pc_L1ZDC_A_C_VTE50`, `HLT_noalg_L1TE600p0ETA49`,
`HLT_noalg_L1TE50_VTE600p0ETA49`), containers `data23_hi.periodAllYear2.physics_{PC,CC}...repro35_v01`.
Pileup removal: **in-time** via a cut on the ZDC–FCal ΣE_T distribution;
**out-of-time** via the ZDC PreSampleAmp (sum over 4 ZDC modules per side),
Gaussian-fit the in-time peak and cut everything > **5σ** from it.

### Glauber-MC parameters (§3, §5.1) — [method-we-use]
- Generator: **PHOBOS GlauberMC v3.2**, **1,000,000** Pb+Pb events at **5.36 TeV**.
- **σ_NN (inelastic nucleon-nucleon):** nominal **68.1 mb** (§3, p.5).
  **⚠ Internal inconsistency in the draft:** the systematics §5.1 (p.16) states the
  varied central value as **67.1 ± 0.5 mb**. The nominal-vs-systematics base values
  differ (68.1 vs 67.1) — flag if this number is quoted; check the approved version.
- **Woods-Saxon (R, a) = (6.62 fm, 0.546 fm)**; systematic ±(0.06, 0.01) fm.
- **Nucleon hard core d_min = 0.4 fm**; systematic ± 0.2 fm.
- **Two-Component Model (TCM):** N_elem = x·N_coll + ½(1−x)·N_part. Single-element
  ΣE_T response fit with Gaussian (chosen) or di-Gamma; event response is an
  N_elem-fold convolution.
  - Fit (Gaussian, chosen): **μ = 11.1 GeV, σ = 5.9 GeV**; (di-Gamma: Θ = 1.95, k = 5.7 GeV).
  - **x = 0.085 ± 0.005** measured; **set to 0.09** to match the previous (5.02 TeV) analysis.
- **Working point:** Glauber matched for **ΣE_T > 30 GeV**, estimated to contain
  **89% ± 1%** of the total Glauber-like cross section. Centrality defined up to
  **85%** (events sorted by ΣE_T, divided into 89 equal-fraction bins).
- Tool: `CentralityFittingTool` / Centrality Analysis Tool (ref [4], CERN GitLab).

### σ_PbPb (total hadronic cross-section) — NOT in this paper
This note gives the **89%±1% Glauber-like efficiency above 30 GeV** and the
Glauber σ_NN, but does **NOT** state the total Pb+Pb hadronic cross-section as a
number in barns. **Our code's σ_PbPb = 7.8 b (roadmap Q2.2) is NOT sourced here**
and remains OPEN — do not attribute 7.8 b to this paper. (Run 2 used 7.66 b at
5.02 TeV; 7.8 b at 5.36 TeV needs its own citable reference.)

### FCal ΣE_T centile cuts (Table 1) — [method-we-use]
Centile cut → FCal ΣE_T threshold. **Values are in TeV** (e.g. the 85% cut =
0.0388 TeV ≈ 38.8 GeV, just above the 30 GeV working point at 89%). More central
= higher ΣE_T. Selected anchor rows (full 1% table in PDF p.10–11):

| Centile | FCal ΣE_T [TeV] | | Centile | FCal ΣE_T [TeV] |
|---|---|---|---|---|
| 85% | 0.0388482 | | 30% | 1.39178 |
| 80% | 0.063208  | | 20% | 2.12188 |
| 75% | 0.0976472 | | 15% | 2.59464 |
| 70% | 0.144347  | | 10% | 3.15972 |
| 60% | 0.288534  | | 5%  | 3.84498 |
| 50% | 0.52171   | | 1%  | 4.51272 |
| 40% | 0.876324  | | 0.9% | 4.53563 |

### ⟨T_AA⟩ per centrality (Table 7) — [method-we-use]  **the number we use**
Units: **mb⁻¹**. (⟨T_AA⟩ = ⟨N_coll⟩ / σ_NN.) Columns: ⟨T_AA⟩, σ(T_AA), σ/T_AA.

**5% intervals (0–80%):**
| Bin | ⟨T_AA⟩ | σ | | Bin | ⟨T_AA⟩ | σ |
|---|---|---|---|---|---|---|
| 0–5%  | **26.1428**¹ | 0.476 | | 40–45% | 3.13726 | 0.126 |
| 5–10% | **20.3241**  | 0.566 | | 45–50% | 2.16697 | 0.139 |
| 10–15%| 15.8262 | 0.652 | | 50–55% | 1.48624 | 0.119 |
| 15–20%| 12.2742 | 0.717 | | 55–60% | 1.00239 | 0.126 |
| 20–25%| 9.62447 | 0.542 | | 60–65% | 0.65588 | 0.111 |
| 25–30%| 7.39028 | 0.343 | | 65–70% | 0.429498| 0.095 |
| 30–35%| 5.56462 | 0.336 | | 70–75% | 0.28466 | 0.078 |
| 35–40%| 4.22415 | 0.257 | | 75–80% | 0.170772| 0.037 |

¹ The PDF text extracts "0–5%" as "326.1428" — an OCR/figure-overlap artifact;
the true value is **26.1428** (matches code; see below).

**10% intervals (0–80%):**
| Bin | ⟨T_AA⟩ | σ | | Bin | ⟨T_AA⟩ | σ |
|---|---|---|---|---|---|---|
| 0–10% | 23.2335 | 0.512 | | 40–50% | 2.65212 | 0.110 |
| 10–20%| **14.0502** | 0.684 | | 50–60% | 1.24431 | 0.121 |
| 20–30%| **8.50738** | 0.441 | | 60–70% | 0.542689| 0.103 |
| 30–40%| 4.89439 | 0.296 | | 70–80% | 0.227716| 0.058 |

**20% intervals (0–100%):** 0–20% 18.6418; 20–40% 6.70088; 40–60% 1.94822;
60–80% 0.385203; 80–100% 0.0560088. **0–1% (1% bins)** range 28.89 (0–1%) → 18.29 (9–10%).

### ⟨N_coll⟩ (Table 5) / ⟨N_part⟩ (Table 6) — selected (full tables in PDF p.18–21)
| Bin | ⟨N_coll⟩ | ⟨N_part⟩ | | Bin | ⟨N_coll⟩ | ⟨N_part⟩ |
|---|---|---|---|---|---|---|
| 0–5%  | 1780.33 | 382.104 | | 30–40% | 333.308 | 129.133 |
| 5–10% | 1384.07 | 330.199 | | 40–50% | 180.609 | 86.364 |
| 10–20%| 956.819 | 258.783 | | 50–60% | 84.7378 | 51.9988 |
| 20–30%| 579.352 | 185.197 | | 60–70% | 36.9571 | 29.1222 |
| | | | | 70–80% | 15.5075 | 15.1449 |

### Systematics (§5)
Two sources, added in quadrature: **geometric** (six alternative Glauber sets for
±1σ up/down of σ_NN, R, a, d_min — see params above) and **efficiency** (vary the
89% efficiency by ±1%, i.e. scale the ΣE_T interval by 90/89 or 88/89). σ/T_AA
ranges from ~1.6% (most central) to ~27% (peripheral 80–100%); see the σ columns.

## Code cross-check — DOES IT MATCH? **YES (all 6 bins).**
Code hardcodes `T_AA = {26.1428, 20.3241, 14.0502, 8.5074, 3.7733, 0.6716}` mb⁻¹
(roadmap Q2.3). These are the **2023 values from Table 7 of this paper**, for the
analysis's 6 centrality classes **0–5, 5–10, 10–20, 20–30, 30–50, 50–80%**:
- 0–5% = 26.1428 ✓ (Table 7, 5% bin)
- 5–10% = 20.3241 ✓
- 10–20% = 14.0502 ✓ (10% bin)
- 20–30% = 8.5074 ≈ 8.50738 ✓
- **30–50% = 3.7733** = mean of Table-7 10% bins (4.89439 + 2.65212)/2 = 3.77326 ✓
- **50–80% = 0.6716** = mean of (1.24431 + 0.542689 + 0.227716)/3 = 0.67157 ✓

So the two wider classes (30–50%, 50–80%) are **simple unweighted averages** of the
paper's 10% bins. (Roadmap 2026-06-15: authoritative 2023 T_AA now stored in
`IntNotes/data/centrality/TaaValues2023.txt`; this paper is their origin.)
The analogous ⟨N_coll⟩/⟨N_part⟩ for those classes would be 30–50%: N_coll 256.96,
N_part 107.75; 50–80%: N_coll 45.73, N_part 32.09 (same averaging).

## References worth future reading   (§6; ≤3)
- **[3] D. Perepelitsa et al., *Centrality determination in √s_NN = 5.02 TeV Pb+Pb
  in 2015*, CERN tech. rep. (CDS 2212936)** — SUPPORTIVE. The direct Run-2
  predecessor whose method (and x = 0.09) this note inherits; gives 5.02 TeV
  centrality/T_AA for comparison. **(Already on disk as `centrality 2015.pdf`.)**
- **[2] Loizides, Kamin, d'Enterria, *Improved MC Glauber predictions...*, PRC 97
  (2018) 054910 [Erratum PRC 99 019901], arXiv:1710.07098** — SUPPORTIVE. Origin
  of the Woods-Saxon / σ_NN / hard-core parameter set used here; useful for
  justifying/varying Glauber inputs in systematics.
- **[1] Loizides, Nagle, Steinberg, *Improved PHOBOS Glauber MC*, SoftwareX 1–2
  (2015) 13** — SUPPORTIVE. The Glauber-MC code itself; cite for the centrality
  tool. (Lower priority — software reference.)

## Related KB docs   (knowledge graph)
- [[atlas_centrality_2015]] — the 5.02 TeV / 2015 predecessor (ref [3]); same FCal-ΣE_T + Glauber method, x=0.09 inherited from there
- [[glauber_modeling]] — Glauber-model fundamentals (N_part/N_coll/T_AA, Woods-Saxon, optical vs MC) underlying this calibration
- [[run2_hf_muon_raa]] — uses the analogous FCal-Glauber centrality + ⟨T_AA⟩ R_AA normalization we adopt
- [[hi_big_picture]] — Glauber/centrality big-picture framing (its "⟨T_AA⟩, σ_PbPb at 5.36 TeV come from elsewhere" pointer resolves to THIS doc)
</content>
</invoke>
