# HIJING Overlay: Truth Pair Structure and HF-Ancestry Limitations

Test sample: Pythia 5.36 TeV pp hQCD pTH8–14 overlaid with HIJING b < 5 fm (ultra-central Pb+Pb).

---

## 1. Truth Pair Structure

### What `truth_mupair_*` actually stores

The overlay NTUP `truth_mupair_bar1/bar2` barcodes are ALL < 200,000 — i.e., **Pythia-only pairs**. The skimmer mupair loop (`TrigRates.cxx:2300`) filters for `id==13 && status==1` from `TruthParticles`. HIJING hard-scatter particles in the overlay `TruthParticles` container are not stable muons (they are quarks, gluons, and hadrons), so only Pythia stable muons form pairs.

Empirically verified: across 200 overlay NTUP events, every `truth_mupair_bar1` and `truth_mupair_bar2` value is < 200,000.

### Why there are so many pairs per event

Overlay event 0 has 78,226 truth particles with bc < 200,000 (vs 486 in PP fullsim). The difference is Geant4 shower secondaries: in the overlay AOD, all Geant4 products (δ-rays, shower electrons, pions, kaons decaying to muons) are stored as truth particles with sequential Pythia barcodes. These include many soft muons from light-hadron decays — they form most of the mupairs.

### Parent group distribution (100-event test, after `truth_parents` fix)

| Group | Fraction |
|-------|---------|
| `s_light` (4) | 68% |
| `resonance_decay` (0) | 29% |
| `direct_b` (1) | 0.5% |
| `direct_c` (3) | 1.8% |

Compare PP fullsim: 51% `direct_b`, 0% `s_light`.

The 68% `s_light` comes from soft Geant4 muons (pion/kaon decays) that dominate the truth-level pair pool. The 0.5% `direct_b` are the genuine Pythia hard-scatter B-decay muons.

---

## 2. `from_same_b` and HF Ancestry Limitations

### Observed counts

| Sample | Events | `from_same_b` pairs |
|--------|--------|---------------------|
| PP fullsim | 10k×6 kn | 93,400 (29.7%) |
| Overlay (test) | 100 | ≈ 0 |

### Root causes

**Cause 1: Geant4 soft muons dominate the pair pool.**  
Most pairs in the overlay flat tree are Geant4 secondary muons from pion/kaon decays. Their parent group is `s_light`, not `direct_b`. Even for the genuine Pythia B-decay muons, they are a tiny fraction of all pairs.

**Cause 2: b-quark absent from overlay truth chain.**  
For Pythia B-decay muons, the ancestry trace shows:
```
μ → B-meson → 3000208 (beam, bc=1)
```
The b-quark is missing (Geant4/AOD truth thinning in overlay suppresses parton-shower intermediates). `FindHeavyQuarks` returns −1 → `skip_event_origin_analysis = true` → `from_same_b` is never set, even for pairs where both muons come from b-decays.

**Note on `truth_parents` reading**: `SetBranchAddress("truth_parents", &ptr)` requires a compiled CollectionProxy for `vector<vector<int>>`. This proxy is only available if `PythiaAnalysisClasses_h.so` (ACLiC pre-compiled) is newer than the source (ROOT auto-loads it), OR via `gInterpreter->GenerateDictionary("vector<vector<int>>","vector")`. Without it, `SetBranchAddress` silently returns −3 and ancestry trace fails. Fixed in `PythiaTruthExtras.c::InitInputExtra` (Apr 2026). The original overlay flat tree (Apr 23) was generated when the header was newer than the .so, so this bug was active.

### Path forward

`from_same_b` efficiency measurements in overlay require:

1. **Apply truth-level pT cut** (≥ 4 GeV) to select hard-scatter Pythia muons before pair formation, eliminating the Geant4 soft-muon contamination.

2. **Replace b-quark matching with B-hadron matching**: modify `HFMuonPairAnalysis` to allow `from_same_b` when both muons trace to the same `m1_eldest_bhadron_barcode` (the B-meson) even without finding the b-quark above it. Currently `skip_event_origin_analysis=true` prevents this.

Both changes together should restore the ~29% from_same_b rate seen in PP fullsim.
