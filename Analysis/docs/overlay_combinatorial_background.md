# HIJING Overlay: Truth Pair Structure and HF-Ancestry Fixes

Test sample: Pythia 5.36 TeV pp hQCD pTH8–14 overlaid with HIJING b < 5 fm (ultra-central Pb+Pb).

---

## 1. Truth Pair Structure

### What `truth_mupair_*` actually stores

The overlay NTUP `truth_mupair_bar1/bar2` barcodes are ALL < 200,000 — i.e., **Pythia-only pairs**. The skimmer mupair loop (`TrigRates.cxx:2300`) filters for `id==13 && status==1` from `TruthParticles`. HIJING hard-scatter particles in the overlay `TruthParticles` container are not stable muons, so only Pythia stable muons form pairs.

### Why there are so many truth particles with Pythia-range barcodes

Overlay event 0 has ~78k truth particles with bc < 200,000 (vs ~486 in PP fullsim). The difference is Geant4 shower secondaries: in the overlay AOD, all Geant4 products (δ-rays, shower electrons, pions, kaons decaying to muons) are stored as truth particles with sequential Pythia barcodes. These include many soft muons from light-hadron decays.

---

## 2. Root Causes and Fixes (resolved Apr 2026)

Three bugs were identified and fixed in `PythiaTruthExtras.c`.

### Bug 1: `truth_parents` branch not read (CollectionProxy missing)

`SetBranchAddress("truth_parents", &ptr)` requires a compiled CollectionProxy for `vector<vector<int>>`. Without it, the call silently returns −3 and ancestry tracing fails entirely (all pairs get `s_light`).

**Fix**: `gInterpreter->GenerateDictionary("vector<vector<int>>","vector")` added to `InitInputExtra`. Also `#include "TInterpreter.h"` added to `PythiaTruthExtras.h`.

### Bug 2: Barcode collision in `GetParticleIndex` cache (critical)

The overlay truth table has ~600–1200 barcode collisions per event: the full `truth_barcode` vector contains the original Pythia muon at index `i1` and a Geant4 secondary (rho⁻, π⁺, φ, etc.) at index `i2 > i1`, both with the same barcode value. The cache was built with `bc_map[bc] = i` (last-writer-wins), so lookups always returned the Geant4 secondary instead of the muon.

With wrong bc_map entries:
- Ancestry tracing started from the wrong particle (non-muon Geant4 secondary)
- Parent groups resolved to `s_light` for nearly all pairs → 68% `s_light`
- `from_same_b` never set → 0%

**Fix**: Changed `barcode_to_index_cache[bc] = i` → `barcode_to_index_cache.emplace(bc, i)`. The `emplace` call is a no-op if the key already exists, so the first occurrence (the Pythia muon) is always preserved.

### Bug 3: `m_eldest_bhadron_barcode` stores beam barcode (not B-meson barcode)

In overlay AODs, the b-quark is absent from the truth chain (AOD truth thinning). The ancestry chain ends at:
```
μ → B-meson → proton (id=2212, bc=1)
```
rather than `μ → B-meson → b-quark → ...` as in PP fullsim.

The B-hadron tracing loop (after bc_map fix) correctly traces into the B-meson, then steps to the parent: a proton (id=2212) at bc=1. At that point, `m_eldest_bhadron_barcode` was set to bc=1 (the **proton's** barcode) for ALL B-meson chains, causing a false `from_same_b` when any two B-meson muons were paired.

**Fix**: When exiting the B-hadron tracing loop, check if the particle above the eldest B-hadron is a b-quark (`abs(first_hadron_id) == 5`). If so (standard MC), store the b-quark barcode as before. If not (overlay: proton/beam), store `last_bhadron_bc` — the eldest B-hadron's own barcode. This correctly identifies genuinely same-B pairs.

---

## 3. Parent Group Distribution After Fixes

500-event test (kn=0, pTH8–14), sign2 (OS pairs):

| Group | Count | Fraction |
|-------|-------|---------|
| `s_light` (0) | 626 / 4096 | 15% |
| `direct_b` (1) | 2132 / 4096 | 52% |
| `b_to_c` (2) | 640 / 4096 | 16% |
| `direct_c` (3) | 624 / 4096 | 15% |
| `from_same_b` | 1228 / 4096 | **30%** |

PP fullsim baseline (same kn=0, 200-event test): `direct_b` 51%, `from_same_b` 30%. Overlay matches closely.

The residual 15% `s_light` comes from Geant4 soft muons (pion/kaon decays) in the overlay truth record. These pairs can be removed by applying the truth-level pT ≥ 4 GeV cut at pair formation (implemented via `PassCuts_PythiaCore` — already present but acts on the individual muon level; the Geant4 muons that pass pT > 4 GeV are unavoidable and constitute genuine combinatorial background).

---

## 4. b-quark Absence and `skip_event_origin_analysis`

For B-meson chains, `FindHeavyQuarks` returns −1 (b-quark not in truth record) → `skip_event_origin_analysis = true`. The "Previous HQ barcode is -1" print is suppressed when `abs(first_hadron_id) != 5` (overlay case). The `from_same_b` check in `HFMuonPairAnalysis` (line 1246) runs **before** the `skip_event_origin_analysis` guard (line 1278), so same-B matching is unaffected by the skip flag.

Pairs with `skip_event_origin_analysis = true` have `pair_origin_analysis_skipped = true` and their `muon_pair_origin_category` is not set. These should be treated as unclassified for origin-level analysis, but their `from_same_b` and `m1/m2_parent_group` fields are valid.
