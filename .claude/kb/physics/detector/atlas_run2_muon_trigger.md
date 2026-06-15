# ATLAS Run 2 Muon Trigger Performance

**Source:** *Performance of the ATLAS muon triggers in Run 2*, The ATLAS
Collaboration, JINST **15** (2020) P09015 — CERN-EP-2020-031.
**arXiv:** 2004.13447
**PDF:** `./2004.13447.pdf` (PRIMARY)
**Classification:** PRIMARY
**Added:** 2026-06-14

## Relevance to this analysis   (REQUIRED)

This is the authoritative reference for the ATLAS muon **trigger** (distinct from
the muon *reconstruction/ID* paper, summarized in
[[ATLAS_Run2_muon_reconstruction]]). Our analysis triggers on muons (PbPb:
single-mu4; pp24: 2mu4) and measures the **trigger efficiency ourselves with a
custom data-driven procedure** (analysis_overview §4b; roadmap steps 8/10), so
this paper is a methodological and numerical **reference**, not a source of
inputs.

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Dimuon trigger efficiency **factorization** ε_dimuon = ε_leg1·ε_leg2·ρ_μμ (§12.4, Eq. p.35) | Directly mirrors our pp24 2mu4 pair efficiency ε_trig^pair = ε₁^nc·ε₂^nc (overview §4b; roadmap step 10; IntNote §10) | [method-we-use] |
| **ρ_ΔR(ΔR, \|y_μμ\|)** close-by-RoI dimuon correction, error-function turn-on → 1 for ΔR≳0.3 (§12.4, Fig. 31) | Conceptual reference for our ΔR trigger-correlation correction ε_dR (currently dummy ≡ 1; roadmap step 9; IntNote §9) | [method-we-use] |
| Single-muon **HLT_mu4 / HLT_mu8** efficiency vs p_T and η, pp + Pb+Pb (§10.2, §11; Figs. 18, 23) | Reference values for our PbPb single-mu4 ε^nc and pp 2mu4 single-leg ε (roadmap steps 8/10; IntNote §8/§10) | [method-we-use] (reference, NOT input) |
| **Tag-and-probe** trigger-efficiency method (J/ψ low-p_T, Z intermediate-p_T) (§8.3) | Reference for the official method vs our custom data-driven one; methodology phrasing | [method-we-use] / [background-for-writing] |
| Trigger **menu / HLT chain naming** (2mu4, mu4, mu8, L1_2MU4 …) (§6, Tables 1–2) | Datasets/trigger section; helps interpret HLT names (IntNote §2/§3; roadmap Q2.4) | [background-for-writing] |
| Overall single-muon trigger ≈68% (barrel)/85% (endcap); HI ≈ pp (Conclusion) | Intro/datasets context on muon-trigger performance | [background-for-writing] |

## Scope & condition-difference warnings

**System / energy / period:** Run 2 (2015–2018). pp at √s = 13 TeV (139 fb⁻¹);
low-pile-up pp at 5.02 TeV (25 pb⁻¹); Pb+Pb at √s_NN = 5.02 TeV (0.49 + 1.42
nb⁻¹); p+Pb at 8.16 TeV. Efficiencies from Z→µµ (intermediate p_T), J/ψ→µµ
(low p_T), W/tt̄ (high p_T).

**Warnings (acknowledge; do NOT quantify size or direction):**
- **Run 2 → Run 3.** Our analysis is Run 3 (pp24 and PbPb 23/24/25) at 5.36 TeV.
  The trigger hardware changed for Run 3 (notably the **New Small Wheel** in the
  endcap, expected to *moderately improve* trigger performance). Run 2 efficiency
  numbers here are reference, not Run 3 values. Do not assume the magnitude of
  the change.
- **Our method ≠ ATLAS official method.** We measure a custom single-muon
  **no-correlation efficiency ε^nc(p_T, q·η)** in our own data and build the pair
  efficiency from it (PbPb inclusion–exclusion ε₁+ε₂−ε₁ε₂; pp 2mu4 product
  ε₁·ε₂). ATLAS here uses Z/J/ψ tag-and-probe and reports scale factors. **Their
  absolute efficiencies/SFs are a cross-check reference, not inputs to our
  correction.**
- pp vs PbPb: higher occupancy in central Pb+Pb; the paper finds HI trigger
  efficiency comparable to pp and **no significant central/peripheral
  difference** for HLT_mu4 (§11), but at Run 2 energies/conditions.

## Content summary

### Trigger system & menu (§5–6)  [background-for-writing]
Two levels: hardware **L1** (six programmable p_T thresholds L1_MU4/6/10/11/15/20;
"11" = 10 GeV threshold with a 3-station RPC coincidence) and software **HLT**
(fast reconstruction → precision stage, RoI-based or full-scan FS). L1 geometric
coverage ≈99% endcap, ≈80% barrel (gaps at η≈0, feet, support structures). HLT
muon: stand-alone (SA, MDT-only) → combined with ID track.

**pp primary chains (Table 1, unprescaled 2016–2018):**
- Single: `HLT_mu26_ivarmedium` (seed L1_MU20, isolated), `HLT_mu50`,
  `HLT_mu60_0eta105_msonly`.
- Multi: `HLT_2mu14` (L1_2MU10), `HLT_mu22_mu8noL1`, `HLT_3mu4(6)`,
  `HLT_3mu6_msonly`, `HLT_mu20_2mu4noL1`.

**HI / low-pile-up chains (Table 2):**
- `HLT_mu8` (Pb+Pb single, seed L1_MU6); `HLT_2mu3` (Pb+Pb dimuon, L1_2MU4);
  `HLT_2mu4` (p+Pb & low-pile-up pp dimuon, L1_2MU4); `HLT_mu14` (low-pile-up pp,
  L1_MU10); `HLT_mu15_L1MU6` (p+Pb high-p_T).
  (Note: chain *names* are Run 2; our Run 3 exact HLT names need user
  confirmation — roadmap Q2.4.)

### Trigger efficiency measurement: tag-and-probe (§8.3)  [method-we-use]
ε = N_match / N_probe, measured **w.r.t. offline-reconstructed muons**. Probe is
"triggered" if within a cone of an object that fired the trigger (ΔR<0.01 for
J/ψ/low-p_T, ΔR<0.1 for Z). Tag = Tight-quality muon that fired the tag trigger.
- **Low p_T (≲10 GeV): J/ψ→µµ**, mass 2.7–3.5 GeV; yields from an extended
  unbinned ML fit (Gaussian + Crystal Ball signal, exponential background;
  simultaneous fit of triggered/not-triggered). To avoid tag–probe correlation,
  reject pairs with ΔR(µ⁺,µ⁻) < 0.2 and impose ΔR(L1 RoI, µ) matching.
- **Intermediate p_T (~10–100 GeV): Z→µµ**, |m−m_Z|<10 GeV, tag p_T>28 GeV +
  Loose iso. **High p_T (≳100 GeV): tt̄ / W+jets** (E_T^miss-triggered tag).
- **Scale factor SF = ε_data / ε_MC** corrects simulation (method biases cancel).
- Systematics (Z): pile-up dependence, tag–probe correlation
  (Δφ(tag,probe)<π−0.1 cut), background (±5 GeV mass window), probe charge/IP/
  isolation/p_T. (J/ψ: statistical only.)

### L1 single-muon efficiency (§9)  [background-for-writing]
Plateau (Z→µµ): barrel L1_MU4 ≈75%, L1_MU20 ≈70% (lower due to 2-/3-station
requirement); endcap L1_MU4 ≈97%, L1_MU20 ≈90% (endcap higher — no barrel
coverage gaps). Barrel turn-on steeper (better p_T resolution).

### pp single-muon HLT efficiency (§10)  [method-we-use: reference]
Primary `mu26_ivarmedium`/`mu50` (Medium muons): **barrel ≈68%, endcap ≈85%** in
26<p_T<100 GeV. Barrel deficit inherited from L1. SFs ≈0.87–0.99 (Table 3); small
p_T dependence; data/MC differ mainly in barrel from RPC modelling.

**HLT_mu4 (Fig. 18) — most relevant to us:** Tight-quality muons, 2015–2016 pp.
**Barrel plateau ≈70%, endcap ≈90%.** Data efficiency **~10% lower than MC in the
barrel**; plateau SF p_T-dependence small. (`HLT_mu14` ≈ similar shape; FS
`HLT_mu8noL1` ≈ 1 because L1 inefficiency is bypassed.)

### Multi-muon (dimuon) trigger efficiency & ΔR correction (§10.2, §12.4)  [method-we-use]
Dimuon efficiency **factorizes** (paper Eq., p.35):

  ε_dimuon = ε_muon(p_T1, q₁·η₁) · ε_muon(p_T2, q₂·η₂) · ρ_μμ(ΔR_μμ, |y_μμ|)

with single-leg ε_muon from J/ψ tag-and-probe (§10.2), and the dimuon
correction further split:

  ρ_μμ(ΔR, |y_μμ|) = ρ_μ(|y_μμ|) · ρ_ΔR(ΔR, |y_μμ|)

- **ρ_μ (asymptotic):** corrects the opposite-charge + dimuon-vertex requirement
  of BLS algorithms. Values **0.983 (barrel), 0.984 (overlap), 0.979 (endcap)**.
- **ρ_ΔR (close-by correction):** inefficiency of resolving two overlapping muon
  RoIs at L1 for closely spaced muons; measured vs ΔR_μμ in three dimuon-rapidity
  regions (barrel |y_μμ|≤1.0, overlap 1.0–1.2, endcap 1.2–2.3); **fitted with an
  error function that approaches unity for ΔR≳0.3** (Fig. 31). This is the direct
  conceptual analog of our pair-level ε_dR ΔR trigger-correlation correction
  (currently a dummy ≡ 1 — roadmap step 9; overview §4b).

> **Connection to our pipeline.** The factorization "product of single-leg
> efficiencies × a ΔR-dependent dimuon correction" is exactly the structure of
> our **pp24 2mu4** correction (ε₁^nc·ε₂^nc, with ε_dR foreseen). Our PbPb
> single-mu4 case instead uses inclusion–exclusion (a pair is triggered if
> *either* leg fires). Do not infer the size of our ε_dR from these Run 2 ρ_ΔR
> curves — different system, energy, and trigger.

### HI single-muon efficiency (§11)  [background-for-writing / reference]
Pb+Pb 5.02 TeV (tag-and-probe, J/ψ & Z): **HLT_mu4** barrel ≈70%, endcap ≈90%,
**no significant central vs peripheral difference**; **HLT_mu8** barrel ≈70%,
endcap ≈85%. p+Pb 8.16 TeV: HLT_mu4 barrel ≈80%, endcap ≈95%. Overall HI trigger
performance **agrees with pp** (paper's headline cross-check).

### BLS (B-physics & Light States) dimuon triggers (§12)  [background-for-writing]
Low-p_T multi-muon chains with **invariant-mass + vertex** requirements (L1Topo
m_μμ and ΔR cuts, e.g. 2<m_μμ<9 GeV for J/ψ/B; Kalman vertex fit at HLT). L1Topo
reduces L1 rate ≈×4 at ~12% HLT-efficiency cost. Context for why low-mass dimuon
triggers exist and how close-by-muon inefficiencies arise (relevant since our
signal is a low-mass, small-ΔR pair). Our analysis does **not** use BLS
mass/vertex trigger requirements — we use plain single-mu4 / 2mu4.

## References worth future reading   (§6)

**None worth adding.** Reference scan: the muon reconstruction/ID reference
([28] *Muon reconstruction performance of the ATLAS detector in pp…*) is already
covered by [[ATLAS_Run2_muon_reconstruction]]. The trigger-menu notes
([39] ATL-DAQ-PUB-2019-001 *Trigger Menu in 2018*; [42–44] 2015–2017 menus) give
exact **Run 2** per-year HLT names/prescales, but our analysis needs **Run 3**
menus (roadmap Q2.4 is to be answered by the user, not by Run 2 menu docs), so
they add no analysis-relevant info. Remaining refs are generic MC-generator/PDF,
detector-TDR, and luminosity papers not specific to this analysis.

## Related KB docs
- [[ATLAS_Run2_muon_reconstruction]] — companion: muon *reco/ID/isolation*
  efficiency (this doc is the *trigger*). Shared method = tag-and-probe (applied
  there to reco/ID, here to the trigger). See its §1 for CB/SA/IO muon types
  referenced above.
- `analysis/run2_dimuon_note.md` — our Run 2 dimuon note (the analysis these
  trigger numbers would feed); trigger-efficiency section there.
