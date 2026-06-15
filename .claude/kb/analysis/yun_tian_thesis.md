# Yun Tian PhD Thesis — Open Heavy Flavor Production in HIC at the LHC

**Source:** Yun Tian, *Open Heavy Flavor Production in Relativistic Heavy Ion
Collisions at LHC*, PhD thesis, Columbia University (2018).
**Report no. / DOI:** CERN-THESIS-2018-078 (defended/dated 27-04-2018).
**PDF:** `./Yun_Tian_thesis-open-HF.pdf` (PRIMARY, committed)
**Classification:** PRIMARY
**Added:** 2026-06-15

> **What this thesis is.** A single document covering FOUR ATLAS HF analyses by the
> same Columbia/Colorado lineage as OUR analysis: (Ch5) HF-muon **R_AA** at
> 2.76 TeV, (Ch6) HF-muon **v_n** at 2.76 TeV — Ch5+Ch6 = the published
> **arXiv:1805.05220** (PRC 98 (2018) 044905); (Ch7) **HF dimuon azimuthal
> correlations** at 5.02 TeV — the *direct predecessor* of OUR Run-2 dimuon
> back-to-back analysis ([[run2_dimuon_note]] / [[run2_dimuon_backtoback_paper]]);
> (Ch8) **γγ→μμ in non-UPC Pb+Pb** at 5.02 TeV = the published
> **arXiv:1806.08708** (PRL 121 (2018) 212301). It carries far more
> implementation detail than any of those papers. **Two distinct values: (1)
> Run-1 method depth for HF-muon R_AA/v_n and the Δp/p fit; (2) a writing
> template for the PI's own Columbia open-HF thesis.**

## Relevance to this analysis   (specific)

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| **Δp/p ("Eloss") momentum-balance template fit** in its *original Run-1 form* (Eq. 5.1, 5.5): the predecessor of the fake-muon purity fit we adapt | Δp/p purity framework — `analysis_overview` §6, roadmap step 16/`task_07`, IntNote §16; mechanism in [[muon_source_template_fits]] | [method-we-use] |
| **HF-muon R_AA** measurement chain (per-event yield / T_AA / pp σ; differential in p_T, centrality) at full method depth | R_AA observable — overview §3, roadmap step 7/`task_06`, IntNote §15 (our numbers are 5.36 TeV; differ) | [method-we-use] (method) / [background-for-writing] (numbers) |
| **Event-Plane + Scalar-Product v_n extraction** (FCal q-vector, resolution, per-(φ−Ψ) Δp/p template fit) for v2/v3/v4 | Flow machinery *if* dimuon v_n is added (overview §3, roadmap §Q2.9 — open scope); 2.76-TeV companion to [[run2_hf_muon_vn]] | [method-we-use] (conditional) |
| **HF dimuon azimuthal-correlation analysis at 5.02 TeV** (Ch7): SS vs OS pairs, combinatorial bkg, broadening vs centrality | The Run-1 predecessor of OUR dimuon analysis; method lineage for [[run2_dimuon_note]] | [background-for-writing] |
| **MC-truth (not tag-and-probe) reco/trigger efficiency** via sigmoid fits, centrality-by-centrality | Contrast for our reco/trigger-eff write-up (our reco-eff is a 3-D PAIR eff) — IntNote §8/§12 | [background-for-writing] |
| **Full Columbia open-HF thesis structure** (chapter flow: QCD/QGP/HF background → detector → ZDC perf → R_AA → v_n → dimuon → γγ → conclusion) | **Template for writing the PI's own PhD thesis** (depth of methodology exposition, how QGP/HF physics is introduced) | [background-for-writing → thesis] |
| **γγ→μμ non-UPC** first measurement (Ch8 = 1806.08708): acoplanarity broadening | γγ is a background we *veto* ([[atlas_gammagamma_mumu_pbpb]]); Intro/thesis context | [background-for-writing] |
| Physics motivation: HF energy loss, dead-cone, R_AA+v_n joint constraint | Intro / thesis ch.1 — defer to [[open_hf_production]] | [background-for-writing] |

## Scope & condition-difference warnings

- **Ch5/Ch6: Run-1, √s_NN = 2.76 TeV** (2011 Pb+Pb, **0.14 nb⁻¹**; pp 2.76 TeV,
  prescale-corrected **≈570 nb⁻¹** = raw 4.1 pb⁻¹ / EF_mu4 prescale 7.037).
  **Single inclusive HF muon**, 4 < p_T < 14 GeV, **|η| < 1** (nominal; |η|<2 for
  some v_n). **No charm/bottom separation** (inclusive HF only — there is **no d0
  fit** here, unlike the later 5.02 TeV papers).
- **Ch7/Ch8: Run-2, √s_NN = 5.02 TeV** (2015 Pb+Pb **0.49 nb⁻¹**; Ch7 pp ref
  **26 pb⁻¹**).
- OUR analysis is **Run 3, 5.36 TeV** (pp24 + PbPb 23/24/25), a **single-b dimuon
  pair** (b→μ + b→c→μ, low mass, small ΔR). ACKNOWLEDGE every energy / era /
  detector / occupancy / trigger-menu difference; **DO NOT estimate its size or
  direction** (GUIDE §5). The T_AA values below are 2.76 TeV — do **not** reuse
  for Run 3 (our values are 5.36 TeV, [[atlas_centrality_2023]]).
- **Reco efficiency here is single-muon and MC-truth-based** (reconstructed/truth,
  sigmoid fit) — **not** tag-and-probe and **not** a pair efficiency. OUR reco
  correction is a 3-D **pair** efficiency ε_reco(pair p_T, pair η, ΔR), not ε₁·ε₂
  (overview §4b; memory `project_pair_reco_effcy_definition`). Do not import.

## Content summary

### Ch5 — HF-muon R_AA at 2.76 TeV  [method-we-use]

**Muon selection (§5.1.3):** Combined muons only; pixel hits ≥1, B-layer ≥1,
SCT ≥7, (pixel+SCT) holes ≤1, |d0_PV|≤5 mm, |z0_PV|≤5 mm, p_ID ≥3 GeV,
p_MS ≥0.1 GeV; plus the **|p_T^MS| > 1.2 GeV** cut (§5.2.2) which removes muons
where the Pb+Pb trigger efficiency < 50% (the trigger sculpts the background at
low p_T^MS — a data/MC mismatch traced to the trigger).

**Discriminant (§5.2.1):** the **momentum balance** ("Eloss", = ∆p/p_ID in
later notation):
> **∆p_loss/p_ID = (p_ID − p_MS − p_param(p_MS,η,φ)) / p_ID**  (Eq. 5.1)
where p_param is the parametrized calorimeter energy loss (preferred to measured
calo energy since muons are non-isolated). Real (HF/prompt) muons peak at ≈0;
π/K-decay-in-flight + hadronic-shower + mis-reco background is broad and shifted
positive. A second discriminant, **scattering-angle significance S** (Eq. 5.2–5.4,
max over scattering centers), was studied but **NOT used** — poor discriminating
power and unexplained centrality dependence. *This is the Run-1 origin of the
Δp/p fit; mechanism in [[muon_source_template_fits]].*

**Template fit (§5.2.3):** dN/d(∆p/p_ID) = N_μ [ f_sig·dP_sig + (1−f_sig)·dP_bkg ]
(Eq. 5.5), 35 bins, χ²/N_DOF with **MINUIT/MINOS** (Eq. 5.6, variances include MC
template weights, Eq. 5.7–5.12). 9 p_T bins (4–14 GeV) × 5 centrality bins
(0-10/10-20/20-30/30-40/40-60%). **Signal templates centrality-dependent;
background templates centrality-integrated** (bkg shown centrality-independent +
limited stats). Templates from **Pythia6 jetjet** (AUET2B, CTEQ6L1): signal from
**muon-filtered** J1–J3 (p_T^μ>3.5 GeV filter), background from core J1–J5; pp uses
JZ1–JZ4 slices. *(Run-1 sample choice; overview §6 notes the Run-3 π/K-enriched MC
equivalent is not yet identified.)* f_sig rises with p_T and toward central
collisions (HF less suppressed than π/K).

**Efficiencies (§5.1.6, §5.1.5):** MC-truth reco efficiency = reco/truth muons,
**sigmoid fit A/(1+e^{α−βt})**, fit centrality-by-centrality. Pb+Pb plateau
**0.768** (≈77%) for p_T≳6 GeV; pp plateau **0.806** (≈81%) from ~5.5 GeV
(Tables 2,3). Trigger EF_mu4_MSonly_L1TE50 (unprescaled), sigmoid fit, **no
centrality dependence**, from min-bias (triggered/all muons).

**Normalization & R_AA (§5.3.2–§5.3.3):** pp differential cross-section
d²σ/dp_T dη from N_corr/(L·Δp_T·Δη); Pb+Pb per-event yield; R_AA formed as the
Pb+Pb per-event yield over ⟨T_AA⟩·σ_pp. T_AA (Table 6, **2.76 TeV** — do not
reuse):

| Centrality | ⟨T_AA⟩ [mb⁻¹] |
|---|---|
| 0–10%  | 23.45 ± 0.37 |
| 10–20% | 14.43 ± 0.30 |
| 20–30% | 8.73 ± 0.26 |
| 30–40% | 5.04 ± 0.22 |
| 40–50% | 2.7 ± 0.17 |
| 50–60% | 1.33 ± 0.12 |
| 60–70% | 0.59 ± 0.07 |

**R_AA systematics (Table 7):** sources = muon-selection cuts (drop pixel/SCT cuts),
**p_T^MS cut** (vary 1.2→0.5/1.5), **π/K background**, **fitting**, **efficiency**.
Largest systematic is **π/K background** (up to ≈11–15% at low p_T); muon cuts
≈3–6%, p_T^MS ≈2–7%, efficiency ≈1.6–2.4%; statistical dominates at 4–4.5 GeV
(≈20–24%). pp luminosity uncertainty **3.1%**.

**Result:** R_AA roughly **independent of p_T**, < 1; **≈0.35 in 0–10% central**
(Abstract). Suppression increases peripheral→central. Compared with theory and
other experiments (§5.3.4).

### Ch6 — HF-muon v_n at 2.76 TeV  [method-we-use (conditional)]

Two methods (broader than the 5.02-TeV [[run2_hf_muon_vn]], which uses EP only):
- **Event-Plane (EP):** FCal q-vector → EP angle Ψ_n; raw v_n from a per-bin
  **Δp/p template fit done in each (p_T, centrality, φ−Ψ_n) bin** to extract the
  HF signal fraction, then v_n^raw from a cosine fit, corrected by the **FCal EP
  resolution Res{nΨ_n}** (data-driven, vs centrality).
- **Scalar-Product (SP):** weights each prompt muon by the q-vector magnitude;
  considered **superior** (robust to poor resolution: v_n^EP → ⟨v_n⟩ for perfect
  resolution, → √⟨v_n²⟩ for poor). Both give **v2, v3, v4** vs p_T and centrality.
Template covariances across φ−Ψ bins handled by **pseudo-experiment** Poisson
resampling. Systematics dominated by the template-fit/background terms (Tables
8–10). *Defer the shared EP machinery to [[run2_hf_muon_vn]].*

### Ch7 — HF dimuon azimuthal correlations at 5.02 TeV  [background-for-writing]

2015 Pb+Pb **0.49 nb⁻¹** + pp ref **26 pb⁻¹**. Intermediate muons 4 < p_T < 8 GeV,
0 < |η| < 2, **Δp/p < 0.15** cut. Azimuthal correlations measured for **same-sign
and opposite-sign** muon pairs; combinatorial background treated explicitly.
Pb+Pb correlations **broader than pp**; **central more broadened than peripheral**
— qualitatively consistent with muon attenuation in the medium. **This is the
direct Run-1 predecessor of OUR dimuon analysis** ([[run2_dimuon_note]],
[[run2_dimuon_backtoback_paper]]); the later Run-2 work extends it (2015+2018,
per-pair efficiency, Δp/p pair-significance purity).

### Ch8 — γγ→μμ in non-UPC Pb+Pb at 5.02 TeV  [background-for-writing]

2015 Pb+Pb **0.49 nb⁻¹**. γγ→μμ identified by small **acoplanarity Aco**
(tight Aco < 0.015 signal region; background d-template from data with Aco>0.02,
asymmetry A>0.15). FCal-ΣE_T centrality. MC-truth efficiencies. **Peripheral:
strong back-to-back correlation** (as in UPC); **central: significant broadening**,
qualitatively consistent with muon attenuation / EM-field effects in the medium.
**This chapter is the published arXiv:1806.08708 (PRL 121 (2018) 212301)** — the
*first* such measurement. *Mechanism/contrast in [[atlas_gammagamma_mumu_pbpb]].*

### γγ→μμ: 1806.08708 vs 2206.12594 (the requested difference)

Both are ATLAS γγ→μμ in **non-UPC** Pb+Pb at 5.02 TeV (NOT the same measurement):
- **arXiv:1806.08708, PRL 121 (2018) 212301** — the thesis's Ch8. **First
  observation/measurement**, 2015 data only (**0.49 nb⁻¹**). Establishes the
  **centrality-dependent acoplanarity broadening** (back-to-back in peripheral →
  broadened in central) as a possible EM probe of the hot medium.
- **arXiv:2206.12594, PRC 107 (2023) 054907** — already in KB as
  [[atlas_gammagamma_mumu_pbpb]]. **Repeats and supersedes** 1806.08708 with **≈4×
  luminosity** (2015 + 2018, **≈1.94 nb⁻¹**), adds **differential cross-sections**
  (unfolded acoplanarity α and k⊥), and **magnetic-field / |Δy| / event-plane
  tests** → finds **no significant magnetic-field effect** on the muons.
This is corroborated by the existing 2206.12594 KB doc, which records that paper
"repeats and supersedes [1806.08708] with 4× luminosity." (Thesis Abstract + Ch8;
[[atlas_gammagamma_mumu_pbpb]] "References / Not added".)

## References worth future reading   (§6; ≤3)

**None new worth adding.** The thesis's own published companions are already
tracked: **arXiv:1805.05220** (PRC 98 (2018) 044905, the 2.76-TeV HF-muon
R_AA+v_n paper = Ch5+Ch6) is already in `FUTURE_READ.md` via
[[run2_hf_muon_raa]]/[[run2_hf_muon_vn]] — and this thesis now supplies its method
depth, leaving only the citeable published version pending. **arXiv:1806.08708**
(γγ, Ch8) is superseded by [[atlas_gammagamma_mumu_pbpb]] (already in KB) and not
worth a separate summary. The remaining bibliography is Run-1-era HF/QGP material
(dead-cone, energy-loss, ALICE/CMS 2.76 TeV) already covered by
[[open_hf_production]] and [[hi_big_picture]].

## Related KB docs   (knowledge graph)

- [[run2_hf_muon_raa]] — 5.02-TeV successor of Ch5; this thesis is the 2.76-TeV
  predecessor with deeper R_AA method exposition.
- [[run2_hf_muon_vn]] — 5.02-TeV successor of Ch6 (EP only); thesis adds the SP
  method and v2/v3/v4 at 2.76 TeV.
- [[muon_source_template_fits]] — concept hub for the Δp/p fit; Ch5 is its Run-1
  origin (momentum balance "Eloss"). No d0 fit here (inclusive HF, no c/b split).
- [[run2_dimuon_note]], [[run2_dimuon_backtoback_paper]] — OUR Run-2 dimuon
  reference; Ch7 (5.02 TeV) is its Run-1 predecessor analysis.
- [[atlas_gammagamma_mumu_pbpb]] — the later γγ→μμ paper (2206.12594) that
  supersedes the thesis's Ch8 (1806.08708).
- [[open_hf_production]] — field-background hub (energy loss, dead-cone, R_AA/v_n).
- [[hi_big_picture]] — QGP/HIC big-picture; thesis Ch2 (QCD/QGP/HF) intro context.
- [[atlas_centrality_2023]] — OUR 5.36-TeV ⟨T_AA⟩ (thesis Table 6 is 2.76 TeV; do not reuse).
