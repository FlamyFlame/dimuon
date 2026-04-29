# HIJING Overlay: Truth Pair Structure and HF-Ancestry Fixes

Test sample: Pythia 5.36 TeV pp hQCD pTH8–14 overlaid with HIJING b < 5 fm (ultra-central Pb+Pb).

---

## 1. Truth Pair Structure

### What `truth_mupair_*` actually stores

The overlay NTUP `truth_mupair_bar1/bar2` barcodes are ALL < 200,000 — i.e., **Pythia-only pairs**. The skimmer mupair loop (`TrigRates.cxx:2300`) filters for `id==13 && status==1` from `TruthParticles`. HIJING hard-scatter particles in the overlay `TruthParticles` container are not stable muons, so only Pythia stable muons form pairs.

### Why there are so many truth particles with Pythia-range barcodes

Overlay event 0 has ~78k truth particles with bc < 200,000 (vs ~486 in PP fullsim). The difference is Geant4 shower secondaries: in the overlay AOD, all Geant4 products (δ-rays, shower electrons, pions, kaons decaying to muons) are stored as truth particles with sequential Pythia barcodes. These include many soft muons from light-hadron decays.

---

## 2. Observed Symptom Before Fixes

| Sample | `from_same_b` rate | `direct_b` parent group |
|--------|--------------------|------------------------|
| PP fullsim | ~30% | ~51% |
| Overlay (before fixes) | **~0%** | **~0.5%** |

Both metrics were near zero in overlay. Root cause: three bugs in `PythiaTruthExtras.c`, described below.

---

## 3. Root Causes and Fixes (resolved Apr 2026)

### Bug 1: `truth_parents` branch not read (CollectionProxy missing)

`SetBranchAddress("truth_parents", &ptr)` requires a compiled CollectionProxy for `vector<vector<int>>`. Without it, the call silently returns −3 and ancestry tracing fails entirely — all pairs get `s_light`, `from_same_b` stays false.

**Fix** (`f748e3f`): `gInterpreter->GenerateDictionary("vector<vector<int>>","vector")` added to `InitInputExtra`. Also `#include "TInterpreter.h"` added to `PythiaTruthExtras.h`.

### Bug 2: Barcode collision in `GetParticleIndex` cache

The overlay truth table has ~600–1200 barcode collisions per event: the `truth_barcode` vector contains the original Pythia muon at index `i1` and a Geant4 secondary (rho⁻, π⁺, φ, etc.) at index `i2 > i1`, both with the same barcode. The cache used `bc_map[bc] = i` (last-writer-wins), so lookups returned the Geant4 secondary rather than the original muon.

Effect: ancestry tracing started from the wrong particle → parent groups resolved to `s_light` for ~68% of pairs → `from_same_b` = 0%.

**Fix** (`dee4050`): Changed to `barcode_to_index_cache.emplace(bc, i)` — `emplace` is a no-op if the key exists, preserving the first (Pythia) occurrence.

### Bug 3: `m_eldest_bhadron_barcode` stores wrong barcode

**Context — what the code was actually doing, and the quark-vs-hadron question:**

The `from_same_b` check in `HFMuonPairAnalysis` compares `m1_eldest_bhadron_barcode == m2_eldest_bhadron_barcode`. The variable name suggests B-hadron level, but examining `SingleMuonAncestorTracing()` before the fix reveals what value was actually stored:

```cpp
// B-hadron tracing loop — steps UP through B-meson chain
while (parent_ids[0] is a B-hadron) {
    first_hadron_barcode = parent_barcode_ABOVE  // UpdateCurParents returns parent's bc
    first_hadron_id      = parent_id_ABOVE
}
// After loop: first_hadron_barcode = barcode of particle ONE LEVEL ABOVE eldest B-hadron
m_eldest_bhadron_barcode = first_hadron_barcode
```

- **PP fullsim (pre-fix)**: The loop exits when the particle above the eldest B-hadron is not a B-hadron. What that particle is depends on which parent appears first in `truth_parents`: often it's a light quark (u, d, s) from the qqbar→B-meson vertex, not the b-quark. So `first_hadron_barcode` and `first_hadron_id` reflect whichever parent was indexed first — not necessarily the b-quark. The old code stored this barcode unconditionally: sometimes b-quark bc, sometimes light-quark bc. Since both muons trace through the same B-meson, the stored barcode is the same for both → `from_same_b` still worked (the value is a coincidental proxy for B-hadron identity, not quark identity).
- **Overlay**: Two classes of B-meson chains exist. (1) **Pythia signal B-mesons** (low bc, e.g. 206, 264): the full ancestry is present — b-quark appears among the B-meson's multiple parents, and `FindHeavyQuarks` finds it normally. These chains are NOT truncated. (2) **HIJING background B-mesons** (bc 3939–73710, measured Apr 2026): the parent chain ends at a beam proton (id=2212, bc=1) instead of a b-quark — HIJING stores less truth ancestry than Pythia. For these, the loop exits when `parent_ids[0]` is the proton → stores `bc=1` (the proton's barcode) for **all** such chains → `m1 == m2 == 1` for all → `from_same_b=true` for all, which is wrong.

  The previous description ("b-quark absent due to AOD truth thinning") was incorrect. The b-quark IS stored for Pythia signal chains; only HIJING-generated B-mesons lack it.

**The user concern** (raised during investigation): "HF tagging is set at the HF hadron level, not the quark level — so the fix shouldn't be needed." This is correct in intent: `from_same_b` *should* work at the B-hadron level. The pre-fix code happened to work for PP fullsim because the stored value (whatever parent came first) was unique per B-meson branch. In overlay it broke because all chains ended at bc=1. The fix makes it explicitly correct at B-hadron level in both samples.

**Fix** (`dee4050`): Track `last_bhadron_bc` (the eldest B-hadron's own barcode) during the tracing loop. After exit:
```cpp
int store_bc = (abs(first_hadron_id) == 5)   // parent is b-quark (rare in practice)
             ? first_hadron_barcode            // store b-quark bc
             : last_bhadron_bc;               // store B-hadron's own bc (default path)
```
Debug output on 100 events (Apr 2026) confirms: **both PP fullsim and overlay now store genuine B-meson PDG IDs** (511, 513, 521, 523, 531, 533 family) at `m_eldest_bhadron_barcode`. The `abs(first_hadron_id)==5` branch is rarely triggered since light quarks typically appear first in `truth_parents` for qqbar→B production. In overlay, `from_same_b=true` pairs (ev=13, ev=99, etc.) share the same B-meson barcode as intended.

---

## 4. Parent Group Distribution After Fixes

500-event test (kn=0, pTH8–14), OS pairs:

| Group | Count | Fraction |
|-------|-------|----------|
| `s_light` | 626 / 4096 | 15% |
| `direct_b` | 2132 / 4096 | 52% |
| `b_to_c` | 640 / 4096 | 16% |
| `direct_c` | 624 / 4096 | 15% |
| `from_same_b` | 1228 / 4096 | **30%** |

PP fullsim baseline (same kn=0): `direct_b` 51%, `from_same_b` 30%. Overlay matches closely.

The residual 15% `s_light` arises from Geant4 soft muons (pion/kaon decays) whose individual pT can pass the ≥ 4 GeV threshold and enter the pair pool. These are genuine combinatorial background in the overlay truth record.

---

## 5. b-quark Absence and `skip_event_origin_analysis`

For HIJING background B-meson chains, `FindHeavyQuarks` returns −1 (b-quark absent from HIJING truth record) → `skip_event_origin_analysis = true`. Measured on 500-event kn=0 test (Apr 2026): **2.2% of B-hadron muons** in overlay hit this path (84/3789); eldest B-hadron barcodes in this subset span 3939–73710. PP fullsim has 0 truncated chains.

The guard `if (abs(first_hadron_id) == 5)` prevents calling `PrintHistory` for these HIJING chains — `first_hadron_id` here is the ID of the particle **above** the eldest B-hadron (the particle that caused the B-hadron tracing loop to exit), not the B-hadron itself. For normal chains it is 5 (b-quark); for HIJING chains it is 2212 (proton). Calling `PrintHistory` on the proton chain would crash: gluon bc=2 has no valid parent in `truth_barcode`.

The warning print "Previous HQ barcode is -1" is suppressed when `abs(first_hadron_id) != 5` (HIJING case, added in `dee4050`).

The `from_same_b` check in `HFMuonPairAnalysis` (line ~1246) runs **before** the `skip_event_origin_analysis` guard (line ~1278), so same-B matching is unaffected by the skip. Pairs with `skip_event_origin_analysis = true` have their `muon_pair_origin_category` unset but `from_same_b` and `m1/m2_parent_group` fields are valid.

---

## 6. Open Questions / Future Investigation

- **Residual 15% `s_light`**: some Geant4 kaon/pion-decay muons survive the ≥ 4 GeV pT cut. Applying an additional truth-level pT cut at pair formation (not just per-muon) may reduce this.
- **`muon_pair_origin_category` for skipped events**: events where the b-quark is absent have no `origin_category` set. Need to decide how to handle these in the origin-level analysis (treat as unclassified or apply separate category).
- **Verification on full dataset**: the 30% `from_same_b` and flavor fractions were measured on a 500-event kn=0 pTH8–14 test. Need to confirm on the full 6-pTH-bin NTP run.
