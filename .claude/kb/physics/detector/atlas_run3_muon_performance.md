# ATLAS Run 3 Muon Performance (Reconstruction + Trigger)

**Primary source (trigger):** *The ATLAS trigger system for LHC Run 3*, The ATLAS
Collaboration, JINST **19** (2024) P06029. **arXiv:** 2401.06630.
**PDF:** `./ATLAS_Run3_trigger_system_JINST19_P06029.pdf` (PRIMARY).

**Source (reconstruction):** *Muon reconstruction performance with the ATLAS
experiment at the LHC using Run-3 pp collision data*, D. Cieri (on behalf of
ATLAS), PIC 2025 proceedings/slides, report **ATL-PHYS-SLIDE-2025-574**.
**URL:** https://cds.cern.ch/record/2947133 ·
**PDF:** `./ATLAS_Run3_muon_reco_PIC2025_SLIDE-2025-574.pdf`.
Underlying ATLAS **public approved plots:** MUON-2023-02 —
https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PLOTS/MUON-2023-02/index.html
(approval refs ATL-COM-PHYS-2024-158 for 2022 data, ATL-COM-PHYS-2025-870 for
2023 data).

**Classification:** PRIMARY (an official published Run-3 performance paper exists
for the trigger; Run-3 reconstruction performance currently exists only as
ATLAS-Preliminary public plots / proceedings — there is **no full Run-3 muon
reconstruction *paper* yet**, the analog of the Run-2 EPJC 81 (2021) 578).
**Added:** 2026-06-15

---

## Relevance to this analysis   (REQUIRED — specific)

Our analysis is **Run 3** (pp24 2mu4; PbPb 23/24/25 single-mu4), so these are the
condition-matched performance references. We measure trigger efficiency ourselves
(custom data-driven ε^nc; overview §4b) and take reco efficiency from MC; these
sources are **references/cross-checks and writing context, not inputs**.

| What in this source | Where it serves OUR work | Use type |
|---|---|---|
| Run-3 muon **reco/ID efficiency** (Z & J/ψ tag-and-probe, 13.6 TeV) at Run-2 level after NSW added to WP definition | Reco-eff context, roadmap step 12 / IntNote §12 (cross-check that Run-3 MC reco eff is sensible) | [method-we-use] (reference, NOT input) |
| **2022 NSW-commissioning** lower endcap data efficiency; recovered from 2023 | Caveat for any pp24/overlay MC vs data reco-eff comparison; datasets/MC section | [background-for-writing] |
| Run-3 **HI (Pb+Pb) trigger menu** thresholds: unprescaled single-muon ~8 GeV; di-muon 4+4 GeV; muon+jet (μ 4 GeV) | Datasets/trigger section, roadmap Q2.4 / IntNote §2; context for our PbPb single-mu4 chain | [background-for-writing] |
| **L1 muon efficiency** (~60% barrel, ~80% endcap for offline Medium, pT>25 GeV) + HLT-rel-to-L1 ~100% | Reference values for our trigger-eff steps 8/10 / IntNote §8/§10 | [method-we-use] (reference, NOT input) |
| **Tag-and-probe** trigger/reco method on Z & J/ψ | Methodology phrasing vs our custom method | [background-for-writing] |
| New Run-3 **inside-out HLT algorithm** for collimated (close-ΔR) di-muons | Context for close-by-muon trigger/reco of our small-ΔR signal pair (steps 9/12) | [background-for-writing] |
| Momentum **scale & smearing / sagitta-bias** calibration scheme (slides) | Detector response / momentum, roadmap step 13 (method demo) | [background-for-writing] |

## Scope & condition-difference warnings

**Trigger paper:** Run-3 pp at **√s = 13.6 TeV**, commissioning year 2022
(performance numbers mostly Data 2022); covers L1+HLT for all signatures incl. a
dedicated **Heavy-Ion (Pb+Pb) trigger-menu** section.
**Reco slides/plots:** Run-3 pp **13.6 TeV**, 2022 (29 fb⁻¹) and 2023 (28 fb⁻¹).

**Warnings (acknowledge; do NOT quantify size or direction — GUIDE §5):**
- **Run 3 differs from Run 2; differences are NOT quantified for our analysis.**
  Hardware changed (notably the **New Small Wheel**, NSW, replacing the endcap
  Small Wheel; new BIS78 RPCs; upgraded L1Muon/MUCTPI/L1Topo). The Run-2 numbers
  in [[ATLAS_Run2_muon_reconstruction]] / [[atlas_run2_muon_trigger]] are not
  Run-3 values. Do not assume the magnitude of any Run-2→Run-3 change for our
  measurement.
- **√s mismatch.** These sources are **13.6 TeV pp**; our analysis is **5.36 TeV**
  pp24 and Pb+Pb. The reco/trigger numbers here are not our operating point;
  acknowledge, do not extrapolate.
- **Our trigger method ≠ ATLAS official.** We build a per-pair efficiency from a
  custom single-muon no-correlation ε^nc(pT, q·η) measured in our own data (PbPb
  inclusion–exclusion ε₁+ε₂−ε₁ε₂; pp 2mu4 product ε₁·ε₂). The tag-and-probe
  absolute efficiencies/SFs here are a cross-check, not inputs.
- **No dedicated Run-3 Pb+Pb muon reconstruction-performance result was found.**
  The Run-3 reco-efficiency plots are **pp 13.6 TeV** only. The trigger paper has
  a Run-3 **HI trigger-menu** section (thresholds/streams) but states HI muon
  triggers "rely on standard reconstruction… and are not discussed further"
  (p.22) — i.e. no separate HI muon-trigger *efficiency* numbers there.

---

## Content summary

> **Dedup note.** The five muon reconstruction types (CB / MS-extrapolated /
> segment-tagged / calo-tagged / inside-out) and the WP scheme (Loose ⊂ Medium ⊂
> Tight, plus High-pT and Low-pT) are unchanged in concept from Run 2 and are
> explained fully in [[ATLAS_Run2_muon_reconstruction]] §1–2. Only Run-3
> specifics are recorded below.

### A. Muon reconstruction & identification (Run 3, reco slides / MUON-2023-02)

[method-we-use: reference]
- **Method:** Z→µµ (probe pT>10 GeV) and J/ψ→µµ (probe pT>5 GeV) **tag-and-probe**;
  ε(X) = N_matched/N_probes, match ΔR<0.05. QCD background subtracted by a
  same-charge template fit in the tag-probe mass spectrum; example background
  fraction ≈ **0.198 ± 0.040%** (Z, two-track probes, bin 5) (slide 7,
  ATL-COM-PHYS-2024-158). Same WP scheme as Run 2: Loose / Medium / Tight /
  High-pT (>300 GeV) / Low-pT.
- **Headline result:** Run-3 identification + reconstruction efficiency is "at the
  level of Run 2, **after the addition of NSW data to the working-point
  definition**" (slide 17, Conclusions). Loose/Medium Z and J/ψ efficiencies sit
  near ~0.95–1.0 across |η|<2.5 (plot axes, slides 8–10); Medium dips around η≈0
  (ATLAS services) due to its tighter n-layers requirement.
- **2022 NSW commissioning caveat:** in 2022 NSW hits were **not used in all
  periods** for WP definition → **lower data efficiency in the endcaps** and
  Data/MC < 1 there; **from 2023 the NSW operates as expected**, bringing Data/MC
  closer to unity (slides 5, 8). ⟨attributed to ATL-COM-PHYS-2024-158 (2022) /
  ATL-COM-PHYS-2025-870 (2023)⟩
- **Momentum calibration (slides 13–16):** opposite-charge (sagitta) bias
  correction applied to data; scale + smearing applied to simulation; sagitta-bias
  maps in (η,φ) slices. Improved alignment in 2023–2024 (NSW + new alignment)
  reduced correction sizes. Data/MC mass agreement near unity after calibration.
  ⟨method only; numbers extract poorly from slides — check PDF if needed⟩

### B. Muon trigger — L1 (trigger paper §3.2)

[method-we-use: reference]
- **L1 threshold nomenclature changed for Run 3:** the quoted threshold is now the
  **50% efficiency point**, not the plateau (as in Run 2). Consequently a Run-3
  **L1 14 GeV single-muon ≈ Run-2 L1 20 GeV**, and a Run-3 **L1 8 GeV di-muon ≈
  Run-2 L1 10 GeV di-muon**, in performance (footnote 2, p.9).
- **Six RPC barrel thresholds** retained (three low-pT 2-station + three high-pT
  3-station); endcap TGC thresholds increased **6 → 15** via the new Endcap Sector
  Logic, with up to 4 quality flags per candidate; MUCTPI now sends up to **32**
  threshold/flag combinations (was 6 in Run 2).
- **L1 efficiency vs offline Medium muons:** about **60% in the barrel** and **80%
  in the endcap** for pT > 25 GeV (stated p.49, citing ref [45]). Barrel L1 maximum
  efficiencies are **~5% lower than at the end of Run 2**, attributed to RPC
  inefficiencies from gas-distribution leaks (§3.2, p.10). L1 charge-ID accuracy is
  ~100% at low pT, degrading with pT (Fig. 4 right).
- **Close-by di-muon recovery (Fig. 5 left):** new RPC Pad/Sector-Logic firmware
  flags two muons in one barrel tower, improving identification of closely spaced
  muons (e.g. boosted J/ψ) — measured vs ΔR_µµ.

### C. Muon trigger — HLT (trigger paper §5.3, §6.2)

[method-we-use: reference]
- **HLT muon reco:** fast (MDT + endcap sTGC/MM, lookup-table pT) → precision
  (offline-like) reconstruction; most chains use **combined** (MS+ID) candidates.
  A **new Run-3 inside-out algorithm** recovers efficiency in **low-pT close-ΔR
  di-muon topologies (e.g. B-meson decays)** by seeding ID-track reconstruction
  from back-extrapolated MS tracks (§5.3, p.41).
- **HLT muon tracking efficiency** > 99% for tracks from the beamline, nearly
  pile-up independent (§5.1.2). MuonIso isolation tracking uses a wider RoI (0.7×0.7,
  φ-restricted to 20 mm in Run 3).
- **Lowest-unprescaled isolated single-muon HLT:** threshold **lowered by 2 GeV to
  pT > 24 GeV** (from Run-2 26 GeV), rate ~200 Hz at L = 1.8×10³⁴ cm⁻²s⁻¹.
  Isolation: scalar ΣpT of tracks in a variable cone (excluding the muon) < **7%**
  of muon pT, cone ΔR = min(10 GeV/pT, 0.3). Support chains: non-isolated
  **pT>50 GeV**; barrel MS-only **pT>60 GeV** (§6.2.1, p.49).
- **Di-muon HLT:** two combined muons **pT>14 GeV each** (unchanged from Run 2),
  rate 24 Hz; plus single-L1-seeded di-/tri-muon chains with sub-leading muons
  pT>8 GeV (di) and (4,4) GeV (tri) reconstructed by the full-scan algorithm.
- **HLT efficiency relative to L1 is close to 100%** in both barrel and endcap
  (tag-and-probe, Z→µµ; §6.2.2, p.49). Measured efficiency is below MC because the
  barrel L1 is simulated with an optimistic lower-bound chamber efficiency; **scale
  factors** from the data/MC difference correct analysis MC.

### D. Heavy-Ion (Pb+Pb) Run-3 trigger menu (trigger paper §4, p.21–22)

[background-for-writing]
- HI menu goal: keep unprescaled pT thresholds as low as possible while limiting
  underlying-event sensitivity. **Unprescaled single-muon trigger pT ≈ 8 GeV**
  (selects semileptonic heavy-flavour decays); **unprescaled di-muon trigger
  pT 4 GeV on both legs** (quarkonium); **muon+jet** trigger with muon pT 4 GeV +
  jet pT 60 GeV (b-jets / soft muons from b-hadron decays).
- Pb+Pb physics streams: `physics_HardProbes`, `physics_UPC`,
  `physics_PC`/`physics_CC` (peripheral/central minimum bias). Calorimeter triggers
  apply a per-layer UE subtraction accounting for azimuthal anisotropy.
- **HI muon triggers use the standard (pp) reconstruction of §5.3/§6.2 and are not
  separately characterized** in the paper (p.22). Exact Run-3 HI HLT chain names /
  prescales for our analysis remain a user-confirmation item (roadmap Q2.4).

### E. NSW impact (cross-source)

[background-for-writing] NSW (new sTGC + MicroMegas, η∈[1.3,2.7], ~100 µm
resolution) replaced the endcap Small Wheel. Reported effects: recovers Run-2-level
endcap reco/ID efficiency once used in the WP definition (reco slides); enables
further L1 endcap rate reduction (trigger paper §3.2). A widely quoted figure is a
**~56% reduction in fake muon trigger rate** from NSW (ATLAS public communications;
⟨not located in either committed PDF — treat as unverified attribution⟩).

---

## References worth future reading   (§6; ≤3)

1. *Studies of the muon momentum calibration and performance of the ATLAS detector
   with pp collisions at √s = 13 TeV*, EPJC **83** (2023) 686, arXiv:2212.07338 —
   **SUPPORTIVE**. Full scale+smearing & sagitta-bias methodology (the slides only
   sketch it). Serves detector response / momentum-scale, roadmap step 13.
   *Caveat: Run-2 13 TeV data.*
2. *The ATLAS Experiment at the CERN LHC: A Description of the Detector
   Configuration for Run 3*, arXiv:2305.16623 — **SUPPORTIVE**. Authoritative NSW /
   BIS78 / muon-system Run-3 hardware description. Serves IntNote datasets/detector
   description (background-for-writing).
3. ATLAS public approved plots **MUON-2023-02** (web) — **PRIMARY-candidate
   webpage**. The full official Run-3 muon reco/ID/isolation efficiency figure set
   (source of the slide numbers). Serves reco-eff cross-check, step 12. *(Partly
   used here via the slides; the full plot set could be summarized when needed.)*

## Related KB docs
- [[ATLAS_Run2_muon_reconstruction]] — Run-2 muon reco/ID/isolation (EPJC 81
  (2021) 578); full muon-type & WP definitions this doc builds on. Run-2 baseline.
- [[atlas_run2_muon_trigger]] — Run-2 muon trigger (JINST 15 (2020) P09015);
  factorized dimuon eff + ρ_ΔR close-by correction; the Run-2 trigger baseline.
- `analysis/run2_dimuon_note.md` — our Run-2 dimuon note (where these efficiencies
  would feed).
