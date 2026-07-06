# ATLAS Muon Quality Working Points (Medium vs Tight)

Reference for the muon-quality selection used in the low-mass dimuon analysis.
Ground truth = MCP `MuonSelectionTool` source in the release used to skim, and
our own bitmask usage in the NTuple processing code.

- **Release:** AthAnalysis **25.2.89** (`x86_64-el9-gcc14-opt`).
- **Tool source:** `.../InstallArea/x86_64-el9-gcc14-opt/src/PhysicsAnalysis/MuonID/MuonSelectorTools/Root/MuonSelectionTool.cxx` (and header `MuonSelectorTools/MuonSelectionTool.h`).
- **Interface / enums:** `.../PhysicsAnalysis/Interfaces/MuonAnalysisInterfaces/MuonAnalysisInterfaces/IMuonSelectionTool.h`.
- **Quality enum:** `.../Event/xAOD/xAODMuon/xAODMuon/versions/MuonEnums.def:113-119` →
  `Tight=0, Medium=1, Loose=2, VeryLoose=3` (lower value = higher purity).

---

## (a) Bitmask → flag mapping (what we AND)

The per-muon integer `muon.quality` is packed **in the skim**, not in Analysis.
Source: `SkimCode/source/HFtrigValidation/src/TrigRates.cxx:1193-1203`, where
`my_quality = m_muonSelection->getQuality(*Muon)`:

| bit (value) | set when | meaning |
|---|---|---|
| `1`   | `muonType()==Combined` | Combined muon (ID track + MS track, global refit) |
| `2`   | `my_quality <= VeryLoose` | passes VeryLoose |
| `4`   | `my_quality <= Loose` | passes Loose |
| `8`   | `my_quality <= Medium` | passes **Medium** (getQuality is Tight or Medium) |
| `16`  | `my_quality <= Tight` | passes **Tight** (getQuality is Tight) |
| `32`  | `passedIDCuts(*Muon)` | ID-track hit requirements (**IDCuts**) |
| `64`  | ME (extrapolated MS) track exists | — |
| `128` | ID track exists | — |
| `256` | `passedMuonCuts(*Muon)` | loose type/author preselection (**MuonCuts**) |

Our analysis selection (verified `Analysis/NTupleProcessingCode/DimuonDataAlgCoreT.c`
~L581-600, `PowhegFullSimExtras.c` L28-40, `PythiaFullSimExtras.c` L36-42):

- **MEDIUM WP** = `quality & (1 | 8 | 32 | 256)` — combined AND getQuality≤Medium AND IDCuts AND MuonCuts.
- **TIGHT WP** = MEDIUM mask with `8`→`16` — i.e. adds the requirement getQuality==Tight.

Because `getQuality()` returns a single enum, "≤Medium" and "==Tight" are decided
by the same function below; the two masks just read different threshold bits.

---

## (b) Medium vs Tight requirements (combined muons)

All of the following are inside `getQuality()` for `muonType()==Combined`
(`MuonSelectionTool.cxx:514-591`). Common preconditions (else → VeryLoose):
author ≠ STACO (L516); `combinedTrackOutBoundsPrecisionHits == 0` (L525); ID and
ME tracks present with `ME cov(q/p) > 0` (L533).

Three discriminating variables are computed:

- **`nprecisionLayers`** — number of precision (MDT/CSC/NSW) layers on the MS track (`fillSummary`, L1638; recomputed at high |η|).
- **`qOverPsignif`** = |q/p(ME) − q/p(ID)| / √(σ²q/p,ID + σ²q/p,ME) — ID/MS charge-over-momentum agreement significance (`qOverPsignificance`, L445-468).
- **`rho'` (ρ′, `rhoPrime`)** = |p_T(ID) − p_T(ME)| / p_T — fractional ID/MS momentum imbalance / momentum-balance variable (L469-478).
- **`reducedChi2`** = primary-track χ² / n.d.o.f. — **normalized track-fit χ²** of the combined fit (L536).

| requirement | MEDIUM (L554-560) | TIGHT (L544-549) |
|---|---|---|
| precision layers | `nprecLayers>1` **OR** (`==1` & `nprecHoleLayers<2` & `|η|<0.1`) | `nprecLayers>1` (no 1-layer exception) |
| q/p significance | `|qOverPsignif| < 7` (flat) [or toroid-off bypass] | `|qOverPsignif| < 7` (flat) **AND** tighter η/p_T-dependent map cut (see below) |
| **normalized track χ²** | **none** | **`reducedChi2 < 8`** |
| ρ′ (momentum balance) | **none** | η/p_T-dependent map cut (see below) |

**Tight adds `passTight(mu, ρ′, qOverPsignif)`** (L1561-1631) — p_T/|η|-dependent cuts
read from calibration histograms in `MuonSelectorTools/.../muonSelection_tightWPHisto.root`:

- kinematic window `p_T ≥ 4 GeV, |η| < 2.5`;
- `4 ≤ p_T < 20 GeV`: `ρ′ < rhoCut(p_T,η)` **and** `qOverPsignif < qOverPCut(p_T,η)` (`tightWP_lowPt_*`);
- `20 ≤ p_T < 100 GeV`: `ρ′ < rhoCut(p_T,η)` (`tightWP_mediumPt_rhoCuts`);
- `100 ≤ p_T < 500 GeV`: `ρ′ < rhoCut` (`tightWP_highPt_rhoCuts`; bin = −1 → no cut);
- `p_T ≥ 500 GeV`: no extra cut.

**Net difference (Medium → Tight):** Tight additionally imposes (i) a **normalized
track-fit χ² < 8**, (ii) **tighter, p_T/η-binned cuts on ρ′ (ID/MS momentum balance)
and q/p significance** than Medium's flat |q/p sig| < 7, and (iii) removes the
single-precision-layer central exception. These are precisely the variables that
detect an ID/MS inconsistency.

### (b′) IDCuts and MuonCuts (bits 32, 256)

- **Combined muon (bit 1):** `muonType()==Combined` — an ID track and an MS track
  matched by a global fit (MuidCB / MuGirl authors). Prerequisite for ρ′ and q/p
  significance to be meaningful.
- **IDCuts (bit 32)** — `passedIDCuts(const TrackParticle&)`, L1478-1508, ID-track hits:
  - pixel hits + dead sensors ≥ 1;
  - SCT hits + dead sensors ≥ 5 (`> 4`);
  - pixel holes + SCT holes < 3;
  - TRT (only `0.1 < |η| ≤ 1.9`): total TRT hits > 5 and outliers < 0.9 × total.
- **MuonCuts (bit 256)** — `passedMuonCuts`, L1459-1476: for combined muons just
  `author != STACO` (a loose type/author preselection; the substantive gates for
  combined muons are IDCuts and the quality WP).

---

## (c) χ² / quality-variable landscape — what is available to cut on

Variables computed by the tool (and hence packable / already in the ntuple stream):

| variable | definition | in Medium? | in Tight? |
|---|---|---|---|
| normalized track χ² (`reducedChi2`) | combined-fit χ²/ndof | no | **yes (<8)** |
| q/p significance | ID/MS charge/momentum pull | flat <7 | flat <7 **+ tight map** |
| ρ′ momentum balance | \|p_T(ID)−p_T(ME)\|/p_T | no | **yes (tight map)** |
| `nprecisionLayers` / holes | MS precision layers/holes | yes | yes (stricter) |
| `combinedTrackOutBoundsPrecisionHits` | out-of-bounds MS hits | veto (=0) | veto (=0) |
| momentum-balance / scattering-curvature / scattering-neighbour significance | used only in the Low-pT MVA WP | no | no |

**Already applied by OUR analysis on top of Medium:** a flat cut
`|dP/P| < 0.12` on every muon (`DimuonDataAlgCoreT.c:599`, threshold
`ParamsSet.h:360 deltaP_overP_thrsh=0.12`). The branch `muon_deltaP_overP` is
defined in the skim as `dP/P = (P_ID − P_ME)/P_ID`
(`SkimCode/.../TrigRates.cxx:1242`) — i.e. **the same physics as ρ′**, the
fractional ID/MS momentum imbalance, but a single flat threshold rather than
Tight's p_T/η-binned map. So we already impose a momentum-balance cut; we do
**not** yet impose any normalized-χ² cut.

---

## (d) Recommendation: would a stricter χ² / q/p cut help remove hadronic & fake muons?

**Physics of the target background.** π/K decay-in-flight muons carry a
**kink** between the ID track (parent hadron) and the MS track (decay muon):
the ID and MS momenta disagree, so these muons show **large ρ′, large q/p
significance, and an inflated combined-fit χ²**. Hadronic punch-through / fakes
give poor ID–MS matching and out-of-bounds MS hits. These are exactly the
variables the Tight WP tightens.

**Assessment.**
1. The **Medium WP does essentially no χ² selection** and only a **loose** flat
   q/p-significance cut (<7). It therefore leaves real headroom: a stricter
   normalized-χ² and/or momentum-balance cut is *not* redundant with Medium and
   would plausibly remove additional decay-in-flight and fake muons in the
   low-mass region.
2. **Tight is the packaged version of exactly this idea** — it adds
   `reducedChi2 < 8` plus p_T/η-binned ρ′ and q/p-significance cuts. Moving
   Medium→Tight, or adding a standalone `reducedChi2 < 8` cut and/or tightening
   our existing `|dP/P|` below 0.12, are all coherent ways to raise purity
   against the hadronic/fake component.
3. **Caveat — most of the low-mass combinatoric background is REAL muons.** The
   dominant low-mass OS background here (single-b, g→QQ̄/flavour-excitation,
   mixed-event combinatoric) is genuine correlated **HF-decay muons**, which are
   well reconstructed and would **not** be removed by a track-quality cut. A
   χ²/ρ′/q/p tightening targets the *fake + decay-in-flight* sub-component only,
   not the physics combinatoric. Its yield impact on the template-fit inputs is
   therefore expected to be modest and must be measured, not assumed.
4. **Downside for signal (real HF-decay muons).** Tight (and any added χ²/ρ′
   cut) costs efficiency, and the loss is **p_T/η-dependent and largest at low
   p_T** — which is our signal region (p_T > 4 GeV, low mass). Our reco
   efficiency is currently a placeholder (see
   `project_pp_reco_eff_placeholder`, `project_pbpb_reco_eff_placeholder`), so
   tightening the WP would require re-deriving the efficiency (and, for full
   Tight, adopting the Tight scale factors) before the change is usable in the
   cross-section / R_AA chain.

**Suggested path (reasoned, not a measured result).** Rather than jump to the
full Tight WP (whose calibration maps and SFs assume standard reco and would
need to be pulled in), study **one transparent, MC-validatable variable at a
time on top of Medium**: (i) add `reducedChi2 < 8`, and/or (ii) scan the
existing flat `|dP/P|` cut (0.12 → e.g. 0.08). For each, compare the low-mass
OS−SS / template yield against the Pythia-fullsim truth-matched **real-muon**
efficiency to confirm the removed events are predominantly fake / decay-in-flight
and that the real-muon efficiency loss is small and modellable. Adopt only if
the purity gain outweighs the efficiency cost and the efficiency change is
re-derived. This is a physics-motivated recommendation grounded in the WP
definitions above, not a claim of measured numbers.
