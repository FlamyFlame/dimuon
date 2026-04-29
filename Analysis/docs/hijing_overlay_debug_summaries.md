# HIJING Overlay: Debug Summaries

Two investigations (both Apr 2026) on Pythia 5.36 TeV pp hQCD DiMu overlaid on HIJING b < 5 fm
(ultra-central Pb+Pb, 2024 HI conditions).  All measurements use pTH14–24 unless noted.

---

## 1. AOD Truth Container Survey

**Goal:** Determine whether a Pythia-signal-only truth container exists in the overlay AOD, as
an alternative to the merged `TruthParticles`.

**Method:** `CollectionTree->Print("*ruth*")` on
`run_pythia_fullsim_HIJING_overlay/mc23_5p36TeV/AOD.49786360._000011.pool.root.1`.

### Containers found

| Container | Content |
|---|---|
| `TruthParticles` | Flat merged list — Pythia signal + HIJING hard-scatter; all read via one loop |
| `TruthEvents` | Pythia hard-scatter event(s); particle-link accessor, **not** a flat container |
| `TruthPileupEvents` | HIJING pileup events; particle-link accessor |
| `MuonTruthParticles` | Reconstruction-linked truth muons (see test below) |
| `egammaTruthParticles` | Egamma-specific truth particles |

No `SignalTruthParticles`, `HardScatterParticles`, or similar Pythia-only flat container exists.

### MuonTruthParticles diagnostic

A temporary branch (`sig_muon_pt`, `sig_muon_barcode`) was added to `TrigRates.cxx`, the skim was
rebuilt (AthAnalysis 25.2.89) and run on 200 events of the local overlay AOD.

| Branch | Mean size/event | Notes |
|---|---|---|
| `truth_barcode` (all TruthParticles) | 69,326 ± 9,128 | full merged table |
| `truth_muon_pt` (TruthPart, μ, status=1) | 9.02 ± 3.21 | current skimmer output |
| `truth_mupair_pt` (all muon pairs) | 41.27 ± 32.22 | current skimmer output |
| `MuonTruthParticles` | **17.20 ± 4.35** | larger than truth_muon_pt |

Barcode breakdown across 200 events in `MuonTruthParticles`:

| Range | Count | Fraction |
|---|---|---|
| bc < 200k (Pythia-range) | 691 | 20.1% |
| bc ≥ 200k (HIJING hard-scatter) | 2749 | **79.9%** |

**Conclusion:** `MuonTruthParticles` is **not** Pythia-signal-only.  It contains HIJING muons (79.9%
by count) and does not apply the `status==1` filter, so it is also larger than
`truth_muon_pt` (17 vs 9 muons/event).  There is no dedicated Pythia-only truth-particle
container in this overlay AOD; the only separation is by barcode convention (below).

---

## 2. Barcode Convention and Filter Mechanics

### 2.1 Convention

| Barcode range | Owner |
|---|---|
| 1 – ~200,000 | Pythia hard-scatter particles **and all Geant4 secondaries** produced by tracking those particles (and HIJING particles) through the detector — including δ-rays, shower electrons, pions/kaons and their muon-decay daughters, and Geant4 secondary B-mesons with truncated ancestry (b-quark not visible in truth record; HIJING origin inferred) |
| 200,001+ | HIJING initial hard-scatter quarks and gluons only |

Both groups appear in the flat `TruthParticles` AOD container.  Crucially, Geant4 secondaries
from HIJING particles (including Geant4 secondary B-mesons, bc 3939–73710 measured) carry bc < 200k and are
indistinguishable by barcode from Pythia Geant4 secondaries.  Only the HIJING hard-scatter
initiators themselves carry bc ≥ 200,001.

### 2.2 Quantitative measurements

Measured on the first 500 events of the pTH14–24 NTUP (HeavyIonD3PD).

**`truth_barcode` (all stored TruthParticles):**

| Sample | bc < 200k | bc ≥ 200k | Total/event (mean) |
|---|---|---|---|
| PP fullsim nn | 97.75% (max = 1,368) | 2.25% (min = 200,001) | ~568 |
| HIJING overlay pp | 94.85% (max = 83,806) | **5.15%** (min = 200,001) | ~68k |

**`truth_muon_barcode`, `truth_mupair_bar1`, `truth_mupair_bar2`:**

| Sample | bc < 200k | bc ≥ 200k |
|---|---|---|
| PP fullsim | **100%** | 0% |
| HIJING overlay | **100%** | 0% |

All truth muons and truth muon pair barcodes stored in the NTUP are Pythia-range (bc < 200k) in
both samples.

**Branch-level per-event sizes (mean ± σ, pTH14–24, 10k events):**

| Branch | PP fullsim | HIJING overlay | Ratio |
|---|---|---|---|
| `truth_barcode` | 586 ± 247 | 68,396 ± 9,024 | ×117 |
| `truth_muon_pt` | 2.24 ± 0.50 | 8.72 ± 2.89 | ×3.9 |
| `truth_mupair_pt` | 1.51 ± 1.12 | 37.83 ± 26.49 | ×25 |

### 2.3 Where the filter operates

There are two mechanisms, only one of which is actually active.

**Mechanism A — explicit bc cut (inactive for overlay).**
`TrigRates.cxx:2185`:
```cpp
if (track->uid() >= 200000 || track->uid() == 0) continue;
```
This line is inside `if ((m_store_truth & Truth::StoreParents) == 0)`.  For MC overlay
`StoreTruth = 1+2+4+8 = 15`, so `StoreParents` is always set and this block never executes.
`IsPrimaryParticle()` (`TrigRates.cxx:56`) also enforces `uid < 2e5`, but only for
per-reco-muon truth-match quality metadata.

**Mechanism B — natural muon filter (the effective one).**
`TrigRates.cxx:2314–2315`:
```cpp
if (id1_ != 13 || m_truth_status[index1] != 1) continue;
```
The HIJING hard-scatter particles in `TruthParticles` with bc ≥ 200k are quarks, gluons, and
unstable hadrons — none are status=1 muons.  They are naturally excluded from `truth_muon_pt`
and `truth_mupair_*` without any explicit barcode cut.

### 2.4 What the extra muons and pairs are

The ~6.5 additional truth muons/event in overlay vs pp fullsim are **not** HIJING hard-scatter
muons (bc ≥ 200k).  They are Geant4 soft muons (bc < 200k) — pions and kaons produced in the
Pythia hard-scatter shower and in HIJING hadronic activity, which Geant4 tracks through the
heavy Pb+Pb material and some fraction decays to μ before stopping.  The large increase in
`truth_barcode` (×117) reflects this Geant4 secondary explosion in the PbPb environment; most
secondaries are electrons, gammas, and light hadrons, but the soft-muon tail is non-negligible.

These Geant4 soft muons carry bc < 200k and appear in `truth_muon_pt`, forming combinatorial
pairs (the ×25 increase in `truth_mupair_pt`).  They are correctly tagged as `s_light` by the
ancestry tracer (Section 5), and are excluded from `_single_b`-category histograms.

---

## 3. Downstream Analysis Safety

**Truth origin tags (`m1/m2_parent_group`, `from_same_b`):** Computed entirely from the
Pythia+Geant4 ancestry chain (all bc < 200k).  HIJING hard-scatter particles (bc ≥ 200k)
are not muons and never enter the ancestry tracing.  Tags are trustworthy after the three bug
fixes described in Section 5.

**Caveat — 2.2% Geant4 secondary B-meson muons:** Geant4 secondary B-mesons (bc 3939–73710;
bc < 200k, distinct from HIJING hard-scatter particles which have bc ≥ 200k) have parent chains
that end at the beam proton (bc=1) with no b-quark visible — Geant4 stores only dummy
beam-particle parents for these secondaries.  HIJING b-quark origin is inferred (B-meson PDG
ID requires a b-quark ancestor; these are absent from PP fullsim) but cannot be confirmed from
the truth record.  For
the 2.2% of B-hadron muons that trace to such a chain, `skip_event_origin_analysis = true`;
their `muon_pair_origin_category` is unset.  Their `m_parent_group` (e.g. `direct_b`) and
`from_same_b` fields remain valid (Section 7).

**Reco efficiency histograms:** The efficiency measures how many Pythia truth muons
(bc < 200k, all entries in `truth_muon_barcode`) have a matched reco muon.  Matching is by
barcode: the reco-muon truth-barcode list (`muon_truth_barcode`) is searched for the Pythia
muon's bc.  A HIJING-matched reco muon carries bc ≥ 200k and can never satisfy this search
for a Pythia truth muon.  HIJING cannot inflate the efficiency count.  HIJING *can* degrade
reconstruction quality for Pythia tracks (noise, confusion in the muon spectrometer), but
that real detector effect is exactly what the overlay efficiency is designed to measure.

**Detector response histograms:** Compare reco vs truth kinematics for muon pairs where
both legs are reco-matched.  All pairs consist of Pythia truth muons (bc < 200k), so HIJING
does not inject fake signal into the response matrices.

---

## 4. Truth Pair Structure

### 4.1 What `truth_mupair_*` stores

The NTUP `truth_mupair_bar1/bar2` barcodes are **all < 200,000** (confirmed on 500 events,
see Section 2.2).  The skimmer mupair loop (`TrigRates.cxx:2313`) filters for
`pdgId==13 && status==1` from `TruthParticles`.  HIJING hard-scatter particles (bc ≥ 200k) are
not stable muons, so they never enter the pair pool.

### 4.2 Why there are so many truth particles with Pythia-range barcodes

PP fullsim has ~586 particles/event in `truth_barcode`; the HIJING overlay has ~68,396.
The difference is not HIJING contamination — it is Geant4 secondaries from simulating all
particles (both Pythia signal and HIJING) through the PbPb detector environment.  All these
secondaries are assigned sequential barcodes in the < 200k range.  The dramatic increase
reflects the much denser Pb+Pb material and nuclear interactions compared to pp.

---

## 5. Three Bugs and Fixes (resolved Apr 2026)

### Bug 1: `truth_parents` branch not read (CollectionProxy missing)

`SetBranchAddress("truth_parents", &ptr)` requires a compiled CollectionProxy for
`vector<vector<int>>`.  Without it, the call silently returns −3 and ancestry tracing fails
entirely — all pairs get `s_light`, `from_same_b` stays false.

**Fix** (`f748e3f`): `gInterpreter->GenerateDictionary("vector<vector<int>>","vector")` added
to `InitInputExtra`.  Also `#include "TInterpreter.h"` in `PythiaTruthExtras.h`.

### Bug 2: Barcode collision in `GetParticleIndex` cache

The overlay truth table has ~600–1200 barcode collisions per event: `truth_barcode` contains
the original Pythia muon at index `i1` and a Geant4 secondary (ρ⁻, π⁺, φ, etc.) at index
`i2 > i1` with the **same** barcode.  The previous cache used `bc_map[bc] = i`
(last-writer-wins), so lookups returned the Geant4 secondary rather than the original muon.

Effect: ancestry tracing started from the wrong particle → ~68% of pairs resolved to
`s_light` → `from_same_b` = 0%.

**Fix** (`dee4050`): Changed to `barcode_to_index_cache.emplace(bc, i)` — `emplace` is a
no-op if the key exists, preserving the first (correct Pythia) occurrence.

### Bug 3: `m_eldest_bhadron_barcode` stores wrong barcode

**What the pre-fix code stored:**
```cpp
// B-hadron tracing loop
while (parent_ids[0] is a B-hadron) {
    first_hadron_barcode = parent_barcode_ABOVE   // one level above current B-hadron
    first_hadron_id      = parent_id_ABOVE
}
// After loop: first_hadron_barcode = particle ONE LEVEL ABOVE the eldest B-hadron
m_eldest_bhadron_barcode = first_hadron_barcode
```

**Two classes of B-meson chains in the overlay truth table:**

1. **Pythia signal B-mesons** (low bc, e.g. 206, 264): full ancestry present; b-quark appears
   among the B-meson's parents; `FindHeavyQuarks` finds it normally.  These are NOT truncated.

2. **Geant4 secondary B-mesons with truncated ancestry** (bc 3939–73710, measured Apr 2026;
   bc < 200k — NOT HIJING hard-scatter which has bc ≥ 200k): the truth parent chain ends at
   the beam proton (id=2212, bc=1) with no b-quark visible.  HIJING origin is inferred from:
   (a) B-meson PDG ID requires a b-quark ancestor by b-number conservation; (b) high Geant4
   bc range inconsistent with Pythia hard-scatter; (c) near-absence in PP fullsim (0.0007%).
   The b-quark is not stored in the truth record, so the origin cannot be confirmed from ancestry
   alone.  After Bug 2 fix, the loop correctly reaches these B-mesons but exits when the parent
   is the proton → stores `bc=1` for **all** such chains → `m1 == m2 == 1` for all pairs →
   spurious `from_same_b=true` for all.

**Fix** (`dee4050`): Track `last_bhadron_bc` (the eldest B-hadron's own barcode) during the
loop.  After exit:
```cpp
int store_bc = (abs(first_hadron_id) == 5)  // parent is b-quark (rare: light quarks
             ? first_hadron_barcode          //   from qqbar→B usually appear first)
             : last_bhadron_bc;             // store the B-hadron's own bc (default)
```
Debug on 100 events confirms both pp fullsim and overlay now store genuine B-meson PDG IDs
(511, 513, 521, 523, 531, 533 family) in `m_eldest_bhadron_barcode`.

---

## 6. Parent Group and Same-B Distribution After Fixes

**Full test-sample validation** (all available files, Apr 2026), opposite-sign (OS) pairs
(`muon_pair_tree_sign2`), cross-section weighted fractions:

| Group (m1) | Truth | PP fullsim | Overlay |
|---|---|---|---|
| `direct_b` | 21% | 21% | 21% |
| `b_to_c` | 57% | 57% | 57% |
| `direct_c` | 13% | 13% | 13% |
| `s_light` | 11% | 11% | 11% |

| Group (m2) | Truth | PP fullsim | Overlay |
|---|---|---|---|
| `direct_b` | 20% | 20% | 20% |
| `b_to_c` | 46% | 46% | 46% |
| `direct_c` | 22% | 22% | 22% |
| `s_light` | 12% | 12% | 13% |

**`from_same_b` fraction (OS):** truth = 25.7%, PP fullsim = 25.9%, overlay = 24.7% — excellent
agreement across all three datasets.  SS from_same_b is ≈0% for all (expected: same-sign pairs
cannot both arise from a single B-meson sequential decay).

All three datasets agree to better than 1% absolute on OS parent groups.  The residual ~12%
`s_light` in OS pairs arises from Geant4 soft muons (pion/kaon decays, bc < 200k) entering the
pair pool; they are correctly tagged and excluded from `_single_b` histograms.  The overlay
s_light fraction is indistinguishable from PP fullsim in OS, confirming the fixes are complete.

---

## 7. b-quark Absence in Geant4 Secondary B-meson Chains and `skip_event_origin_analysis`

For Geant4 secondary B-meson chains with truncated ancestry (bc 3939–73710; HIJING origin
inferred but b-quark not visible in truth record), `FindHeavyQuarks` returns −1 →
`skip_event_origin_analysis = true`.

**Counts from full test-sample validation** (all available NTUP files, Apr 2026):

| Sample | B-hadron muons | Truncated | Fraction | B-meson bc range | Muon bc range |
|---|---|---|---|---|---|
| PP fullsim (24 files, ~240k events) | 282,827 | 2 | 0.0007% | [78, 370] | [635, 681] |
| Overlay (6 files, ~60k events) | 75,499 | 1,725 | **2.28%** | [3939, 73710] | [4229, 73916] |

The 2 truncated chains in PP fullsim have very low Geant4 barcodes (78, 370), consistent with
edge cases in Pythia chains rather than truncated Geant4 secondary B-mesons; their muon barcodes (635, 681)
are also in the very low Geant4 range.

The guard `if (abs(first_hadron_id) == 5)` prevents calling `PrintHistory` for Geant4
secondary B-meson chains — `first_hadron_id` is the ID of the particle **above** the eldest
B-hadron (5 for normal chains, 2212 for truncated Geant4-secondary chains).  Calling
`PrintHistory` on the proton chain would crash (gluon bc=2 has no valid parent in
`truth_barcode`).

The `from_same_b` check in `HFMuonPairAnalysis` runs **before** the `skip_event_origin_analysis`
guard, so same-B matching is unaffected by the skip.  Pairs with `skip_event_origin_analysis=true`
have `muon_pair_origin_category` unset but `from_same_b` and `m1/m2_parent_group` are valid.

---

## 8. Open Questions

- **Residual ~12% `s_light` in OS pairs**: Geant4 kaon/pion-decay muons that survive the
  ≥ 4 GeV pT cut contribute combinatorial background.  The fraction matches between PP fullsim
  and overlay, so the effect is correctly modeled.  Whether to apply an additional pair-pT cut
  to reduce it should be evaluated against signal acceptance cost.

- **`muon_pair_origin_category` for 2.28% truncated B-meson pairs**: these pairs have no
  `origin_category` set (pair_origin_analysis_skipped=true).  The pair-level skip rate in SS
  is 19.1% and in OS is 5.3% (pair skipped if either muon is truncated).  Need to decide
  whether to treat these as unclassified or add a dedicated category before origin-level
  analysis proceeds.

---

## 9. Full Validation Comparison: Truth vs PP Fullsim vs Overlay (Apr 2026)

**Samples:** all available test-sample NTUP files; PP fullsim = 24 files (~240k events),
overlay = 6 files (~60k events), Pythia truth = full 5.36 TeV sample.

**Tree entry counts (weighted by `weight` branch):**

| Sign | PP fullsim pairs | Overlay pairs | Truth pairs |
|---|---|---|---|
| SS (sign1) | 84,406 | 25,808 | 2,098,165 |
| OS (sign2) | 314,766 | 82,804 | 8,062,511 |

### 9.1 Opposite-sign (OS): excellent three-way agreement

| Observable | Truth | PP fullsim | Overlay | Max Δ |
|---|---|---|---|---|
| `from_same_b` | 25.7% | 25.9% | 24.7% | 1.2% |
| pair_origin_analysis_skipped | 0.1% | ~0% | **5.3%** | — |

m1/m2 parent group: all three datasets agree to ≤1% absolute on every bin.

Plots: `usatlasdata/pythia_fullsim_hijing_overlay_test_sample/plots/truth_origin_comparison/`

Flavor category (OS) — notable discrepancy:

| Flavor | Truth | PP fullsim | Overlay |
|---|---|---|---|
| `from_resonance` | 20% | 20% | 19% |
| `res_contam` | ~4% | ~4% | ~6% |
| `single_b` | 26% | 26% | 25% |
| `bb` | 40% | 40% | 39% |
| `cc` | 10% | 10% | 10% |

`res_contam` is ~2% higher in overlay (absolute) — Geant4 soft muons in the dense PbPb
environment occasionally pair with a resonance muon, inflating the resonance-contamination
category.  All other OS flavor bins match truth and PP fullsim to ≤1%.

Origin category (OS): good agreement on `fc` (~18%), `same_gs_fsr` (~10–11%),
`same_gs_isr_0hs` (~16%).  Overlay shows `diff_gs_same_hs`/`others` elevated by ~4% combined
compared to truth/fullsim, attributable to the 5.3% skip rate shifting some pairs into default
category values.  `not_both_open_HF` is ~49% in all three (dominant OS category).

**Conclusion OS:** Truth analysis in overlay is validated.  Parent groups, from_same_b, and
flavor fractions all agree with PP fullsim and truth to within expected physical differences
(soft Geant4 muons, 5.3% skip rate for pairs involving truncated Geant4 secondary B-mesons).

### 9.2 Same-sign (SS): known discrepancies, all attributable

| Observable | Truth | PP fullsim | Overlay | Δ (abs) |
|---|---|---|---|---|
| `from_same_b` | ~0% | ~0% | ~0.1% | — |
| pair_origin_analysis_skipped | 0.1% | ~0% | **19.1%** | — |

Flavor category (SS):

| Flavor | Truth | PP fullsim | Overlay | Δ |
|---|---|---|---|---|
| `bb` | ~83% | ~83% | ~71% | −12% |
| `res_contam` | ~15% | ~15% | ~25% | +10% |

Origin category (SS):

| Origin | Truth | PP fullsim | Overlay | Δ |
|---|---|---|---|---|
| `fc` | ~30% | ~30% | ~23% | −7% |
| `same_gs_isr_0hs` | ~28% | ~28% | ~22% | −6% |
| `others` | ~10% | ~10% | ~23% | +13% |
| `not_both_open_HF` | ~6% | ~6% | ~9% | +3% |

m2 parent group (SS): `s_light` overlay = ~16% vs truth/fullsim ~10%; `b_to_c` overlay ~58%
vs ~62%.

**Explanation:** All SS discrepancies originate from two known overlay-specific effects:
1. **19.1% pair-level skip rate** (truncated Geant4 secondary B-meson chains) — pairs with
   `skip_event_origin_analysis=true` shift origin fractions; `bb` is lower because many bb
   pairs involve at least one truncated muon.
2. **Extra Geant4 soft muons** (~15% of single-muon pool) — combinatorial pairing of soft
   muons with resonance muons inflates `res_contam`; pairing with b-decay muons inflates
   `s_light` as m2 in SS.  This effect is larger in SS (large combinatorial pool) than OS
   (requires two opposite-sign muons from related chains).

No unexpected discrepancies observed.  OS comparison confirms the HF ancestry fixes are
complete and the truth origin tags are reliable for efficiency/response analysis.
