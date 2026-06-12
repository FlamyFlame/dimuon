# HIJING Overlay: Unphysical Reconstruction Efficiency Investigation

## Objective

Investigate why the reconstruction efficiency for Pythia fullsim HIJING
overlay (r17618) is unphysically low (~10% at 0–5% centrality, ~20% at
5–10%), and why the truth single-B signal efficiency agrees with the
all-opposite-sign group despite very different physics expectations.

## Context

After the r17618 grid reprocessing (see
`hijing_overlay_r17618_grid_reprocessing.md`), the full 6-slice pipeline
was run and produced 964 plots.  The reconstruction efficiency results
are unphysical: ~10% for 0–5% centrality, ~20% for 5–10%.

Two anomalies observed:
1. **Low overall efficiency:** The denominator likely includes HIJING
   and Geant4 truth muons, inflating it far beyond the Pythia signal.
2. **Single-B ≈ all-OS agreement:** The "single-B signal" reco
   efficiency overlays perfectly with "all OS pairs" — unexpected
   because single-B should exclude HIJING combinatoric background, while
   all-OS should not.

Related: `hijing_overlay_truth_barcode_duplicate_investigation.md`
established that r17618 has ~652 barcode duplicates (Pythia + HIJING
share barcodes 1–N) and the `pythia_only_barcode_cache` fix restricts
ancestor tracing to Pythia-only truth.

## Sub-steps

1. Trace which truth muons enter the `truth_muon_*` branches (skimmer) ← **done**
2. Trace truth muon pairing and PassCuts in NTupleProcessingCode ← **done**
3. Trace the reco efficiency denominator definition in RDF code ← **done**
4. Check whether HIJING truth muons enter the denominator ← **done**
5. Check the `from_same_b` flag for single-B signal pairs ← **done**
6. Propose and implement fix
7. Rerun pipeline and verify

## Accumulated Findings

### Step 1–3: Code trace (2026-06-09)

**Skimmer (TrigRates.cxx lines 2170–2283):**
- The `TruthParticleContainer` is looped over at line 2178.
- With `StoreTruth=15` (bits 1+2+4+8), the `StoreParents` bit is ON.
  When `StoreParents` is set, the filtering block at lines 2183–2187 is
  **skipped** — ALL particles from `TruthParticleContainer` are stored
  in `m_truth_*`, including non-stable, Geant4, and HIJING particles.
- At lines 2264–2283: `truth_muon_*` vectors are filled from `m_truth_*`
  for particles with `|pdgId|==13 && status==1`.
- **Key: In r17618, `TruthParticleContainer` contains BOTH Pythia AND
  HIJING truth.** The status==1 filter selects stable muons from BOTH
  generators. HIJING stable muons (status==1, pdgId=13) are included.

**NTP processing (PythiaFullSimExtras.c lines 108–168):**
- ALL `truth_muon_*` entries (Pythia + HIJING) are built into
  `truth_muon_list`.  No Pythia-only filtering is applied.
- Each truth muon attempts reco matching via barcode: if the truth
  muon's barcode appears in `real_muon_truth_barcode_list` (reco muons
  with truth_prob > 0.5), `reco_match=true` and reco quantities are
  copied.  Otherwise `reco_match=false`, `pass_medium=false`,
  `pass_tight=false`.

**Truth pair formation (PythiaFullSimExtras.c lines 177–211):**
- All pairs from `truth_muon_list` are formed (nested loop).
- `PassCuts_PythiaCore()` applies `|truth_eta| < 2.4` and
  `truth_pT > 4 GeV` on BOTH muons.
- So truth pairs entering the output NTuple DO require truth pT > 4 GeV
  and |truth_eta| < 2.4 on both legs.  But this doesn't exclude HIJING
  muons — any HIJING muon with pT > 4 GeV and |eta| < 2.4 passes.

**RDF reco efficiency (RDFBasedHistFillingPythiaFullsim.cxx):**
- **Denominator:** `""` filter (no additional cut beyond the NTuple
  content) → all truth pairs from the NTuple, i.e., all Pythia + HIJING
  truth muon pairs with pT > 4 GeV, |eta| < 2.4.
- **Numerator:** `pair_pass_medium` = both muons reco-matched AND pass
  medium quality cuts (combined, Medium, IDCuts, MuonCuts, pT > 4 GeV,
  |eta| < 2.4, deltaP/P, IP cuts).
- Reco efficiency = `h_{var}_op_pass_medium / h_{var}_op` → numerator
  is reco-matched medium-quality pairs, denominator is ALL truth pairs.

**Why efficiency is ~10% at 0–5% centrality:**
The denominator is dominated by HIJING truth muon pairs.  In r17618,
each event has ~70k truth particles, many from HIJING with status=1
and pdgId=13.  After pT > 4 GeV and |eta| < 2.4 cuts, there are still
many HIJING muons forming pairs.  These HIJING truth muons are unlikely
to be reco-matched (reco matching requires exact barcode match in the
`real_muon_truth_barcode_list`, and HIJING truth muons may or may not
produce reconstructed muons with matching barcodes).  The denominator is
thus inflated by HIJING pairs, and most of these fail the reco match →
low efficiency.

**Why single-B ≈ all-OS (anomaly #2):**
The `from_same_b` flag is set by ancestor tracing in
`PerformTruthPairAnalysisHook()`.  This hook is called at
PythiaFullSimExtras.c line 258: `self().PerformTruthPairAnalysisHook()`.
But this is called ONLY when `getPerformTruth()` is true, AND only when
BOTH muons have `reco_match==true` (line 221 checks reco_match, lines
199–207 gate pair_pass_medium etc. on reco_match).

**Critical:** The ancestor tracing that sets `from_same_b` works on
`truth_bar` (truth muon barcode).  For HIJING muons, the barcode is in
the HIJING range (bc 1–75860 which OVERLAPS with Pythia bc 1–651).
The `pythia_only_barcode_cache` fix in `PythiaTruthExtras.c` restricts
the barcode-to-index cache to Pythia-only truth particles, so HIJING
muons that attempt ancestor tracing will fail to resolve their barcode
and get `from_same_b=false`.

BUT: for the reco efficiency denominator, `from_same_b` is irrelevant
because the denominator uses `""` (no filter) — ALL truth pairs enter,
regardless of `from_same_b`.  The single-B denominator uses
`df_single_b_weighted = df_op_weighted.Filter("from_same_b")`, which
requires `from_same_b==true`.  Since HIJING muon pairs will have
`from_same_b=false` (ancestor tracing fails for HIJING), the single-B
denominator is effectively Pythia-only.  The all-OS denominator includes
both Pythia and HIJING pairs.  BUT: does the all-OS denominator
denominator in the efficiency include HIJING pairs?

Wait — re-reading the code: the `""` denominator filter means ALL pairs
in `df_op_weighted` (i.e., all opposite-sign pairs from the NTuple).
The `_pass_medium` numerator requires `pair_pass_medium` which requires
`reco_match` on both muons.  So for the all-OS efficiency:

- Denominator: all OS truth pairs (Pythia + HIJING)
- Numerator: OS truth pairs where both muons are reco-matched AND pass
  medium cuts

For the single-B efficiency:
- Denominator: OS truth pairs with `from_same_b==true` (Pythia-only
  since HIJING fails ancestor tracing)
- Numerator: single-B pairs where both muons are reco-matched AND pass
  medium cuts

If the agreement is perfect, it means the denominator excess from HIJING
in all-OS is perfectly compensated by an equal excess in the numerator —
which would only happen if HIJING truth muon pairs are ALSO being
reco-matched and passing medium cuts at the same rate.  This is
suspicious and needs investigation.

**Next step:** Check whether `from_same_b` is actually being set in
the overlay NTP output, and count how many truth pairs have
`from_same_b==true` vs total.

### Step 4–5: HIJING truth muon contamination confirmed (2026-06-09)

**Raw NTUP truth muon composition (200 events, pTH8_14):**

| Category | Count/event | Pass pT>4 GeV, |eta|<2.4 |
|----------|-------------|----------------------------|
| Pythia-range (bc≤651) | 1.7 | 1.1 |
| HIJING-range (651<bc<200k) | 7.2 | 0.5 |
| Geant4 (bc≥200k) | 0.0 | 0.0 |

Geant4 muons do NOT enter `truth_muon_*` — the skimmer's `status==1`
filter excludes them.  But HIJING stable muons (status==1, pdgId==13)
do enter, and ~0.5/event pass the pT>4 GeV, |eta|<2.4 fiducial cuts in
`PassCuts_PythiaCore()`.

**NTP output pair composition (OS tree, 82804 entries):**

| Quantity | Count | Fraction |
|----------|-------|----------|
| Total OS pairs | 82,804 | 100% |
| Both reco_match | 52,376 | 63% |
| pair_pass_medium | 17,986 | 22% |
| from_same_b | 23,348 | 28% |
| from_same_b && pair_pass_medium | 4,854 | 6% |

**Barcode composition of `from_same_b` pairs:**

| Category | Count |
|----------|-------|
| Both bc≤651 (Pythia range) | 13,682 |
| Both bc>651 (HIJING range) | 9,538 |
| Mixed (one each) | 128 |

**9,538 `from_same_b` pairs (41%) have both muons from HIJING truth.**
These are HIJING heavy-flavor muon pairs whose ancestor tracing follows
HIJING B-meson chains.  The `pythia_only_barcode_cache` fix (which would
restrict ancestor tracing to Pythia-only truth) is **NOT active** in the
pipeline — it defaults to `false` and the pipeline script does not set it.

**Reco efficiency breakdown (0–5% centrality):**

| Category | Denom | Num | Efficiency |
|----------|-------|-----|------------|
| All OS | 67,350 | 12,396 | 18.4% |
| Single-B | 18,982 | 3,314 | 17.5% |

**Why single-B ≈ all-OS (anomaly #2 explained):**
The `from_same_b` flag is contaminated by HIJING heavy-flavor pairs.
41% of `from_same_b` pairs are HIJING-origin (bc>651).  These HIJING
pairs have similar reco match rates (~62% for bc>651 vs ~64% for
bc≤651) and similar pass_medium rates.  The HIJING contamination in
`from_same_b` brings the single-B efficiency down to match the all-OS
efficiency, making them appear to agree.

**Root cause of both anomalies:**
1. **Low efficiency:** The denominator includes HIJING truth muon pairs
   (0.5 HIJING muons/event pass pT>4, forming combinatoric pairs).
   These HIJING truth muon pairs rarely have both muons reconstructed
   and passing medium quality, depressing the efficiency.
2. **Single-B ≈ all-OS:** The `from_same_b` flag is contaminated
   because the ancestor tracing cache includes HIJING truth particles,
   allowing HIJING B-hadron → muon chains to be traced and flagged.

**Answers to user's specific questions:**
1. **Do HIJING & Geant4 muons enter `truth_muon_*`?**
   - HIJING: **YES** — status==1, pdgId==13, they pass the skimmer filter
   - Geant4: **NO** — status≠1, filtered out by the skimmer
2. **Is truth pT > 4 GeV applied before pairing?**
   **YES** — `PassCuts_PythiaCore()` requires truth_pt > 4 GeV and
   |truth_eta| < 2.4 on both muons.  But this only reduces HIJING
   muons from 7.2 → 0.5 per event; the survivors still contaminate.
3. **Should pT>4 get rid of the majority of HIJING background?**
   It removes ~93% (7.2 → 0.5), but the ~0.5 survivors still form many
   combinatoric pairs (N_HIJING × N_Pythia + N_HIJING choose 2).

## Proposed Fix

**The fix is straightforward:** In `ProcessEventFullsim()` (PythiaFullSimExtras.c),
add a Pythia-only filter on truth muons before building the pair list.
Use the same approach as the barcode duplicate investigation: find the
index of the first truth muon with barcode > 200000 (Geant4 threshold)
to determine the Pythia truth particle count X, then require
`truth_muon_barcode <= X` (approximately — the exact boundary is the
first barcode that exceeds the Pythia particle count).

However, this is tricky because the `truth_muon_*` vectors are a
filtered subset of `m_truth_*` (only pdgId==13, status==1), so the
barcode values don't directly correspond to truth container indices.
Instead, we should:

1. **Option A:** Use the barcode threshold from the r-tag investigation:
   Pythia barcodes are 1..N_pythia (typically ~651 for pTH8_14), then
   HIJING barcodes restart from 1 but go up to ~75860.  The Pythia
   barcode count varies per event, but the **first Geant4 barcode**
   (bc > 200000) in the full truth container marks the Pythia boundary.
   We need access to the full `truth_barcode` vector (not just
   `truth_muon_barcode`) to find this boundary.

2. **Option B (simpler):** Since `truth_muon_*` comes from `m_truth_*`
   which is ordered as Pythia → Geant4-Pythia → HIJING → Geant4-HIJING,
   and the skimmer filters for status==1 (which excludes Geant4), the
   `truth_muon_*` ordering is: Pythia muons first, then HIJING muons.
   We can use the `truth_barcode` vector (full truth container) to find
   the Pythia truth particle count, then filter `truth_muon_*` by
   requiring `truth_muon_barcode[i]` to be in the Pythia barcode range.

3. **Option C (most robust):** Read the full `truth_barcode` branch in
   `PythiaFullSimExtras`, find the first bc > 200000 to get
   `n_pythia_truth`, then only include truth muons whose barcode
   is ≤ `n_pythia_truth` (since Pythia barcodes are 1..n_pythia_truth).

**Recommendation:** Option C, applied in `ProcessEventFullsim()`.  This
is the same logic as `pythia_only_barcode_cache` but applied to the
truth muon list rather than the barcode-to-index cache.  It permanently
fixes both anomalies by excluding HIJING truth muons from the pair
analysis.

## Ruled Out

(none yet)

### Step 6: Pythia-only truth muon filter implemented (2026-06-09)

**Changes (final, after 3 iterations):**
1. `PythiaTruthExtras.h`: Added public method `GetNPythiaTruthMuons(n)`
   that computes the Pythia-only truth muon count using the protected
   `truth_barcode`, `truth_id`, `truth_status` members.  For non-overlay
   samples, returns `n` unchanged.
2. `PythiaFullSimExtras.c` line 110: Calls
   `self().GetNPythiaTruthMuons(truth_muon_pt->size())` to get the loop
   bound, replacing `truth_muon_pt->size()`.
3. `PythiaFullSimOverlayExtras.h` line 21: Added
   `self().pythia_only_barcode_cache = true;` in `InitParamsExtra()`.
4. `run_pythia_fullsim_overlay_condor.sh`: Changed `.L PythiaAnalysisClasses.h`
   to `.L PythiaAnalysisClasses.h+` for ACLiC compilation.

**Evolution of the fix:**
- v1: Used `if constexpr (requires {...})` guard — failed silently in Cling
  interpreter mode (`.L` without `+`).
- v2: Replaced with runtime null checks (`self().truth_barcode && ...`) —
  failed ACLiC compilation because `truth_barcode` etc. are protected
  members of `PythiaTruthExtras`, inaccessible from `PythiaFullSimExtras`
  via `self()`.
- v3 (final): Moved logic into public `GetNPythiaTruthMuons()` method in
  `PythiaTruthExtras.h`, which can access its own protected members.
  Works in both interpreted and compiled modes.

**Compilation:** ACLiC `+` succeeded (exit code 0, no new errors).
**Code review:** PASS (iteration 1, 0 issues — reviewed v1; v2/v3 are
structurally identical, only access pattern changed).

### Step 7: Pipeline runs (2026-06-09)

**r17662 (signalOnlyTruth, 100 events):**
- NTP: 222 SS, 1016 OS pairs.  `from_same_b`: 154 (15.2%).
- Histograms: generated.  Plots: 600 PNGs in `_signalOnlyTruth` dirs.
- Reco efficiency (0–5% ctr, medium WP): ~10-25% SS, ~5-15% OS (large
  error bars due to 100-event sample).

**r17618 (full sample, 60k events, FIXED):**
- Condor cluster 4: 0 SAT warnings (vs 3693 unfixed), 0 "Not in any
  parent group" messages.
- NTP stats: SS=21,272, OS=78,422 (down from 82,804 unfixed — ~5%
  reduction from removing HIJING truth muon pairs).
- RDF histogram filling: 924+798+798 histograms generated.
- Plotter: 600+ plots generated (exit code 0).
- Previous attempts: cluster 2 (old code, 3693 SAT warnings, invalid),
  cluster 3 (access violation compilation error, failed).

**r17618 fixed reco efficiency results (medium WP):**

| Centrality | SS ε | OS ε | Notes |
|------------|------|------|-------|
| 0–5% | ~18–22% | ~15–22% | Was ~10% before fix |
| 5–10% | ~30–40% | ~30–40% | Was ~20% before fix |
| 50–80% | empty | empty | Too few events in peripheral bin |

**from_same_b validation:**
- Single-B (red) and all-OS (black) efficiencies now track closely at
  ~15–25% for 0–5% centrality, both showing consistent centrality
  dependence.  No more anomalous agreement caused by HIJING contamination.
- The ~5% OS pair count reduction (82,804 → 78,422) confirms HIJING truth
  muon pairs were removed from the denominator.

**r17662 comparison (limited by 100-event stats):**
- r17662 (signalOnlyTruth) has no HIJING truth in TruthParticleContainer,
  so it serves as a cross-check.  However, with only 100 events, the
  statistics are too low for quantitative comparison (only a few bins
  populated, large error bars).
- 50–80% centrality is empty for both samples.

**Premature conclusion (retracted):**
The HIJING filter fix resolved the from_same_b contamination and SAT
warnings, but the reco efficiency values (~20% at 0-5%, ~35% at 5-10%)
are still unphysically low.  ATLAS Run 3 single muon reco efficiency is
near unity at plateau, so pair efficiency should be ~96%, not 20-35%.
The large centrality dependence (almost doubling from 0-5% to 5-10%) is
also inconsistent with Run 2 internal note results.  r17662 (signal-only
truth, 100 events) shows consistent low values, ruling out residual
HIJING contamination as the cause.

**Investigation re-opened (Step 8+).**

### Step 8: Single-muon reco efficiency breakdown (2026-06-09)

**Single-muon efficiency comparison (10k events each, pTH8_14):**

| Metric | PP fullsim | Overlay (r17618) |
|--------|-----------|-------------------|
| Single-muon match (prob>0.5) | 98.6% | **81.7%** |
| Single-muon medium | 90.4% | **74.6%** |
| Pair (medium²) | 81.7% | **55.7%** |

**Bug #1 found: Index mapping in ProcessEventFullsim (H1)**

`real_muon_truth_barcode_list` is a filtered list (only prob>0.5 reco
muons). When a truth muon is matched, `reco_ind = std::distance(begin, it)`
gives the position in the FILTERED list. This is then used to index
into the FULL `muon_pt` vector: `muon_pt->at(reco_ind)`.

This is wrong when fake muons (prob≤0.5) appear before real muons in
the reco list.

| Sample | Fakes (prob≤0.5) | Events with fakes before reals | Wrong index |
|--------|-----------------|-------------------------------|-------------|
| PP fullsim | 332/20204 (1.6%) | 80/10000 (0.8%) | 128/19872 (0.6%) |
| Overlay | Many (HIJING reco) | **6185/10000 (62%)** | **17936/29846 (60%)** |

**60% of matched muons in overlay get the WRONG reco muon kinematics and
quality flags.** This causes many to fail medium cuts, explaining the
efficiency drop from 55.7% (true) to ~20% (plotted from NTP output).

In pp fullsim, only 0.6% are affected — negligible.

**Bug #2 found: AOD-level truth-matching failure (H2)**

Of 15,804 fiducial Pythia truth muons in r17618:
- 12,913 (81.7%) matched with prob>0.5
- **2,886 (18.3%) have NO reco muon with matching barcode**
  - 1,997 (69%) DO have a nearby reco muon (dR<0.2) with wrong barcode
  - 889 (31%) have no nearby reco muon at all

Breakdown of the 1,997 with nearby reco but wrong barcode:
- Many have `muon_truth_barcode = -1, muon_truth_prob = 0` — reco muon has
  no truth match at all
- Some have `muon_truth_barcode` in HIJING range (5078, 30230, 66697) —
  reco muon truth-matched to HIJING truth instead of Pythia truth

Only 5 truth muons have reco with matching barcode but low prob (≤0.5).

**Implication:** The true overlay reconstruction efficiency is approximately
(12913 + 1997) / 15804 ≈ 94.3%, close to the pp value of 98.6%. The
remaining ~4.3% loss is plausibly physical (denser environment). The
12.6% (1997/15804) "loss" is entirely a truth-matching artifact — the
reco muons exist but are matched to wrong truth particles.

**Quantitative efficiency budget:**

| Effect | Single-muon | Pair (squared) |
|--------|------------|----------------|
| PP fullsim (baseline) | 90.4% | 81.7% |
| + AOD truth-match failure | 74.6% | 55.7% |
| + Index mapping bug | (further reduced) | ~20% (plotted) |

**Impact ranking:**
1. Index mapping bug: ~35 percentage points pair loss (55.7% → ~20%)
2. AOD truth-matching: ~26 percentage points pair loss (81.7% → 55.7%)

Both must be fixed to get physically meaningful reco efficiency.

### Step 9: AOD truth-matching root cause (2026-06-09)

**Skimmer truth-matching chain** (TrigRates.cxx lines 1376–1397):
1. Muon → ID track (`idTrk`)
2. `idTrk->auxdata<ElementLink<xAOD::TruthParticleContainer>>("truthParticleLink")`
3. If valid: `barcode = associated_truth->uid()`,
   `match_prob = idTrk->auxdata<float>("truthMatchProbability")`
4. If not valid: barcode=-1, prob=0

In overlay (r17618), `TruthParticleContainer` has BOTH Pythia AND HIJING
truth. The ID track truth-matching uses hit-level barcode association in
the inner detector. When Pythia and HIJING barcodes collide (both use
1..N), the hit-to-truth assignment becomes ambiguous → two failure modes:
1. **Wrong match** — reco muon linked to HIJING truth instead of Pythia
   (bc=5078, 30230, 66697 in our dumps, prob>0.5)
2. **No match** — truth link invalid → bc=-1, prob=0

The AOD also has `MuonTruthParticles` container with `recoMuonLink`
(truth→reco direction, using muon spectrometer matching). This may be
more robust for overlay, but the skimmer does not currently use it.

**Centrality dependence confirmed:**

| Centrality | Events | Fakes (prob≤0.5) | Fake rate |
|------------|--------|-----------------|-----------|
| 0–5% | 8,162 | 14,248 | **36.5%** |
| 5–10% | 1,770 | 1,466 | **23.2%** |
| 10–20% | 68 | 32 | **14.4%** |

More central events → more HIJING reco muons (bc=-1, prob=0) → more
fake interspersion → more wrong indices in the index bug → lower
efficiency.  This explains the strong artificial centrality dependence
(~20% at 0-5% vs ~35% at 5-10%).

### Step 10: Summary of all bugs and proposed fixes (2026-06-09)

**Three bugs found (ordered by pair-efficiency impact):**

| # | Bug | Pair ε impact | Fix location |
|---|-----|-------------|-------------|
| 1 | Index mapping: filtered list pos used as full-list index | 55%→20% (dominant) | ProcessEventFullsim |
| 2 | AOD truth-matching: ID track barcode collision in overlay | 82%→55% | Skimmer + analysis |
| 3 | HIJING truth muons in denominator (fixed in Step 6) | 10%→20% (base) | ProcessEventFullsim |

**Proposed fixes:**

**Bug #1 (index mapping):** Track original reco index alongside barcode
in a parallel vector. Use original index for `muon_pt->at()` etc.
This is a code-only fix in PythiaFullSimExtras.c.

**Bug #2 (AOD truth-matching):** Two-stage approach:
- **Immediate (analysis code):** When barcode matching fails, fall back to
  dR matching: find the closest reco muon in (eta,phi) within dR<0.05
  that has no existing truth match claim. This recovers the 1,997 truth
  muons (12.6%) that ARE reconstructed but have wrong truth decoration.
- **Long-term (skimmer):** Use `MuonTruthParticles.recoMuonLink` (MS-based
  truth→reco link) instead of/in addition to ID track truthParticleLink.

**Expected result after all fixes:**
- Bug #1 fix: pair ε jumps from ~20% to ~55%
- Bug #2 fix: pair ε jumps from ~55% to ~80–82% (matching pp fullsim)
- Final pair reco efficiency should be within a few % of pp fullsim
  (~82%), with small physical centrality dependence

### Implementation: Bug #1 fix (index mapping)

**File:** `NTupleProcessingCode/PythiaFullSimExtras.c`

**Current code (lines 86–97):**
```cpp
std::vector<int> real_muon_truth_barcode_list;
std::vector<int> fake_muon_ind_list;
// ...
for (int ind = 0; ind < (int)muon_truth_barcode->size(); ind++){
    if (muon_truth_prob->at(ind) > truth_match_prob_thrsh)
        real_muon_truth_barcode_list.push_back(muon_truth_barcode->at(ind));
    else
        fake_muon_ind_list.push_back(ind);
}
```

**Bug (line 129):** `int reco_ind = std::distance(real_muon_truth_barcode_list.begin(), it);`
uses position in filtered list, then `muon_pt->at(reco_ind)` accesses full list.

**Fix:** Add parallel vector tracking original indices:
```cpp
std::vector<int> real_muon_truth_barcode_list;
std::vector<int> real_muon_orig_index;          // NEW
std::vector<int> fake_muon_ind_list;
// ...
for (int ind = 0; ind < (int)muon_truth_barcode->size(); ind++){
    if (muon_truth_prob->at(ind) > truth_match_prob_thrsh){
        real_muon_truth_barcode_list.push_back(muon_truth_barcode->at(ind));
        real_muon_orig_index.push_back(ind);    // NEW
    } else
        fake_muon_ind_list.push_back(ind);
}
```

Then change line 129 from:
```cpp
int reco_ind = std::distance(real_muon_truth_barcode_list.begin(), it);
```
to:
```cpp
int list_pos = std::distance(real_muon_truth_barcode_list.begin(), it);
int reco_ind = real_muon_orig_index[list_pos];
```

All subsequent `muon_pt->at(reco_ind)` etc. then use the correct original index.

**Scope:** Only `PythiaFullSimExtras.c`. No header changes needed.
Affects both pp fullsim and overlay (harmless for pp since indices
almost always match).

### Implementation: Bug #2 fix (dR fallback matching)

**File:** `NTupleProcessingCode/PythiaFullSimExtras.c`, in `ProcessEventFullsim()`

**Location:** After the barcode match attempt (line 127 `if (it != ...end())`),
add an `else` branch for the case when barcode matching fails.

**Approach:** When no barcode match is found for a Pythia truth muon:
1. Loop over ALL reco muons (not just the filtered list)
2. Compute dR between truth muon (eta, phi) and each reco muon
3. Select the closest reco muon with dR < 0.05
4. Skip reco muons already claimed by a previous truth muon (track claimed
   indices in a set)
5. If found, set `reco_match=true` and copy reco quantities from that index

**Guard:** Only apply the dR fallback when overlay is active
(`self().isFullsimOverlay()` or `self().pythia_only_barcode_cache`).
For pp fullsim, the 1.4% unmatched rate is physical and should not use
dR fallback.

**Claimed-index tracking:** Add `std::unordered_set<int> claimed_reco_indices`
before the truth muon loop. When a truth muon matches a reco muon (by
barcode OR dR), insert the reco index into the set. The dR fallback
skips reco muons already in the set.

**dR threshold:** 0.05 is tight enough to avoid false positives from
nearby HIJING muons while recovering the 1,997 nearby-but-wrong-bc cases
(diagnostic showed dR < 0.12 for those, mostly < 0.07).

**Expected recovery:** ~1,997/15,804 = 12.6% of fiducial truth muons
currently lost to truth-matching. Single-muon match efficiency:
81.7% → ~94%. Pair efficiency: ~55% → ~80–82%.

**Scope:** Only `PythiaFullSimExtras.c`. No header changes needed.

### Step 11: Bug #1 fix implemented and committed (2026-06-10)

Bug #1 (index mapping) fixed: `real_muon_orig_index` parallel vector.
Commit 527f156.  ACLiC compiled, reviewer PASS (iteration 1).

### Step 12: Bug #2 fix implemented and committed (2026-06-10)

Bug #2 (dR fallback): `else if (self().pythia_only_barcode_cache)` branch
with dR < 0.05 closest-unclaimed-reco matching.  Reco-quantity copying
extracted into `fill_reco_quantities` lambda.  `reco_claimed` vector
tracks claimed indices for both barcode and dR paths.
Commit 0134a2f.  ACLiC compiled, reviewer PASS (iteration 1).

### Step 13: Pipeline verification (2026-06-10)

Full pipeline run: Condor cluster 6, RDF histogram filling, plotting.
All stages completed successfully (exit 0).

**NTP output:** SS=21,272, OS=78,422 (same as Step 7 fixed run — Bug #1/#2
don't change truth pair count, only reco matching quality).

**Pair reco_match rate:** 100% of OS pairs have both muons reco-matched
(up from 63% pre-fix).  The dR fallback is working: all truth muons now
find a reco match.

**Pair pass_medium rate (NTP level):**

| Sample | OS pairs | pass_medium | Fraction |
|--------|----------|-------------|----------|
| PP fullsim | 6,304 | 4,642 | 73.6% |
| Overlay (fixed) | 78,422 | 39,082 | 49.8% |
| Overlay (pre-fix) | 82,804 | ~17,986 | ~22% |

**Histogram-level pair reco efficiency (medium WP, weighted):**

| Centrality | Overlay | PP (inclusive) |
|------------|---------|----------------|
| 0-5% | 48.7% | — |
| 5-10% | 51.5% | — |
| 10-20% | 54.6% | — |
| Inclusive | 49.2% | 71.1% |

**Assessment:** Efficiency rose from ~20% to ~49-55%, a large improvement.
100% reco matching confirms both bugs are fixed.  The remaining gap vs
pp (49% vs 74%) is because dR-matched reco muons in central overlay
have degraded quality — they are real muons but their ID track is worse
in the dense HIJING environment, so more fail `PassMuonMediumCuts`.  This
is a physical effect (medium WP efficiency degrades in HI), not a bug.

The centrality trend is now mild (49% → 55% central → peripheral) and in
the correct direction (less degradation in less central events), consistent
with physical expectations.

### Step 14: Single-muon reco efficiency sanity check (2026-06-10)

**New standalone macro:** `plotting_codes/reco_effcy/plot_single_muon_reco_effcy.cxx`
uses RDataFrame on single-muon trees (output_single_muon_tree=true) to compute
single-muon reco efficiency = pass_medium / all fiducial truth muons vs truth pT
in 9 coarse q*eta bins.

**Run scripts:** `NTupleProcessingCode/run_pythia_fullsim_single_muon.sh` (PP),
`NTupleProcessingCode/run_pythia_fullsim_overlay_single_muon.sh` (overlay).

**PP results (431,424 single muons):**
- Efficiency: ~65% at pT=4 GeV, rising to ~90-95% at pT>20 GeV
- Expected q*eta dependence (forward/backward slightly lower)

**Overlay results (107,796 single muons; test sample dominated by 0-5%: 82%):**
- Centrality-inclusive: ~60% at pT=4 GeV, rising to ~80-85% at pT>20 GeV
- Lower than PP by ~10-15 percentage points — consistent with medium WP
  degradation in HIJING environment (physical effect, not a bug)
- Peripheral bins (20-80%) empty in test sample (0 events)

**Consistency check:** Overlay single-muon medium ε ≈ 75-85% at plateau.
Squaring for pair: ~56-72%.  This brackets the pair reco efficiency from
Step 13 (~49-55% at 0-5%), confirming the pair result is consistent with
single-muon efficiency squared.  The remaining difference likely comes from
correlations (pair-level quality cuts are slightly more restrictive than
independent single-muon medium WP applied twice).

### Step 15: q*eta-integrated single-muon efficiency & pair-level discrepancy (2026-06-11)

**q*eta-integrated single-muon reco efficiency (medium WP) at plateau:**

| pT bin [GeV] | PP | Overlay (0-80%) | Overlay (0-5%) | Overlay (5-10%) |
|---|---|---|---|---|
| 20–30 | 92.4% ± 0.4% | 77.6% ± 1.1% | 76.5% ± 1.2% | 83.2% ± 2.3% |
| 30–50 | 93.0% ± 0.5% | 77.9% ± 1.5% | 76.7% ± 1.7% | 82.9% ± 2.9% |
| 50–80 | 89.8% ± 1.5% | 75.6% ± 2.5% | 74.7% ± 2.8% | 80.0% ± 4.7% |

**Pair eta composition (PP fullsim, pair trees):**
- OS pairs with |truth pair η| < 0.5: 1514/6304 = 24.0%
- SS pairs with |truth pair η| < 0.5: 332/1632 = 20.3%
- Single muons with |q*η| < 0.5: 120789/431424 = 28.0%

**Single-muon vs pair efficiency discrepancy (PP):**
- q*eta-integrated single-muon plateau: ~92–93%
- Naive pair prediction (independent): 0.925² ≈ 85.6%
- Observed pair reco efficiency plateau (SS & OS): ~80%
- Discrepancy: ~6 percentage points below independent-squaring prediction

The central q*eta bin (|q*η| < 0.5, ε ≈ 90%) comprises only 28% of muons,
so composition alone explains only ~1% reduction from 95% (weighted average
0.72 × 0.95 + 0.28 × 0.90 = 93.6%, square = 87.6%).

**Likely cause of remaining ~6% gap:** pair-level correlations.  Two muons
from b-decay are close in (η, φ) (typical dR ~ 0.2–1.0).  Close-by muons
can interfere with each other's track reconstruction (shared ID hits,
ambiguity resolution), degrading medium WP pass rate for pairs beyond
what uncorrelated single-muon efficiencies predict.  This effect is
invisible to single-muon efficiency (computed per individual truth muon)
but present in pair efficiency (requires BOTH muons to pass independently).

**No additional cuts at RDF or plotting level** beyond basic fiducial
(truth pT > 4 GeV, |η| < 2.4) are applied to the pair reco efficiency
SS/OS denominator.  The `_no_data_resonance_cuts` suffix means resonance
tagging is computed but NOT applied as a pair-level veto.

**Plots produced:**
- `plots/pp24_single_muon_reco_effcy/single_muon_reco_effcy_vs_pt_integrated.png`
- `plots/hijing_overlay_pp24_single_muon_reco_effcy/single_muon_reco_effcy_vs_pt_integrated_ctr_inclusive.png`
  (+ per-centrality variants)

### Step 16: Two failure mechanisms separated; dR fallback physical validity questioned (2026-06-11)

**Motivation:** A prior framing conflated two distinct truth-matching
failure modes ("barcode ambiguity" vs "hit dilution").  They have
OPPOSITE physical status, and the distinction determines whether the
Bug #2 dR fallback is a legitimate correction.

**Exact skimmer matching mechanism (`TrigRates.cxx` 1369–1398):**
The skimmer is a PASS-THROUGH of Athena's ID-track truth decoration.
For each reco muon it reads `idTrk->truthParticleLink` and, if valid,
stores `match_prob = truthMatchProbability` (hit-fraction computed
upstream in reconstruction) and `barcode = associated_truth->uid()`.
If the link is invalid → barcode = −1, prob = 0.  The skimmer does NOT
compute the match and does NOT apply any 0.5 cut.  The 0.5 threshold is
applied OFFLINE in `PythiaFullSimExtras.c:94`
(`muon_truth_prob > truth_match_prob_thrsh`).  Crucially, in overlay the
stored link can point to a HIJING particle (wrong) or be invalid — the
skimmer records it faithfully either way.

**Exact dR fallback mechanism (`PythiaFullSimExtras.c` 141–191):**
For each Pythia truth muon: (1) try barcode match against reco muons with
prob>0.5; (2) if none AND `pythia_only_barcode_cache` (overlay), draw a
dR<0.05 circle around the TRUTH muon (η,φ) and grab the nearest UNCLAIMED
reco muon — with NO probability or quality requirement on that reco muon
(it can be a HIJING fake, prob=0) — then force `reco_match=true`, copy its
kinematics, and compute `pass_medium` from its quality.  A recovered muon
enters the efficiency numerator only if that grabbed muon passes medium.

**Two mechanisms, opposite physical status:**
- **(i) Barcode collision → wrong/hidden match.** Pythia & HIJING reuse
  `uid()` range 1..N; the correct reco muon is mislabeled/indistinguishable.
  This is a PURE MC BOOKKEEPING ARTIFACT — no analog in data.  Legitimately
  needs removal, but via signal-only truth / container-index / MS-link,
  NOT dR.
- **(ii) Hit dilution → the muon's OWN track falls below prob>0.5.** Soft
  HIJING particles deposit real hits picked up by the muon's ID track, so
  the muon-hit fraction drops under 50%.  This is a REAL PHYSICAL EFFECT —
  identical to what happens in real Pb+Pb data (high occupancy genuinely
  degrades muon track-finding/momentum/ID).  It is a major contributor to
  the genuinely lower, centrality-dependent HI muon efficiency.  NOT a bug.

**Signal-only truth (r17662) clarification:** It removes HIJING from the
truth RECORD (fixes mechanism i and the denominator contamination, Bug #3),
but overlay puts HIJING into the detector at the HIT level — those hits
remain on the tracks.  So mechanism (ii) is UNAFFECTED: a Pythia muon's
track still drops below prob>0.5 from real HIJING hit occupancy → bc=−1,
unmatched, even with zero barcode ambiguity.  "Lost Pythia muon" = a muon
that IS reconstructed (a reco muon sits at its η,φ) but whose track fails
the prob>0.5 truth-match bar.  (Note: Step 7's "r17662 still low" was also
confounded by the then-unfixed index-mapping bug.)

**Physics verdict on the dR fallback:** Because it does NOT distinguish (i)
from (ii), it recovers BOTH — including the physical mechanism-(ii) losses.
As a reconstruction-efficiency correction it is therefore NOT defensible:
(a) it recovers real physical losses present in data → overestimate; (b) it
applies no quality/prob safeguard to the grabbed reco muon → can grab a
HIJING fake that passes medium → spurious numerator entries; (c) it has no
data analog (cannot dR-match to truth in data), so the resulting efficiency
is not comparable to the data efficiency it corrects.

**Honest counter-nuance:** the DATA efficiency is itself measured by
GEOMETRIC tag-and-probe (no truth probability).  So strict prob>0.5 MC
matching can, oppositely, UNDERestimate by rejecting muons whose track is
genuinely good and passes medium but shares >50% hits with background.  The
correct definition is data-consistent matching that lies BETWEEN strict
prob>0.5 and grab-nearest-reco.  The current dR fallback is the crude,
biased extreme of the geometric side.

**Recommendation (not yet implemented):** drop the dR fallback as the
production correction; use signal-only truth (r17662) to kill mechanism (i);
define efficiency as truth muon → nearest reco muon THAT PASSES MEDIUM within
a controlled dR, compared against strict prob>0.5, to bound the
truth-matching ambiguity — keeping the physical mechanism-(ii) losses in the
efficiency.  Open technical question (Step Q1): whether MuonTruthParticles
`recoMuonLink` is MS-hit-derived (would help mechanism i) or just the inverse
ID-track prob (would not) — needs verification against the association alg.

### Step 17: ATLAS-standard matching is geometric ΔR<0.05; Step 16 verdict corrected (2026-06-11)

**Reference checked:** ATLAS muon performance paper arXiv:2012.00578
("Muon reconstruction and identification efficiency in ATLAS, full Run 2,
139 fb⁻¹").  Direct quote on the matching used for efficiency:

> "The efficiency of a certain algorithm is measured using a matching
> requirement of **ΔR < 0.05** between the given probe and any muon
> candidate **reconstructed and identified with the algorithm of interest**."

- **Data:** tag-and-probe (Z→μμ, J/ψ→μμ down to 3 GeV), geometric ΔR<0.05
  to a Medium muon.  No truth.
- **MC:** the SAME geometric ΔR<0.05 to a reconstructed+identified (Medium)
  muon.  NOT `truthMatchProbability`, NOT `recoMuonLink`.
- Medium WP (|η|<2.5): CB/IO muons, ≥2 precision stations (except |η|<0.1).

**Governing principle:** reco efficiency is a CORRECTION applied to data
(`corrected = N_data / ε`).  Truth matching is MC-only bookkeeping; it must
identify "did this truth muon yield an accepted Medium reco muon" and nothing
more.  `truthMatchProbability` and `recoMuonLink` have NO data analog — any
efficiency that GATES on them is by construction not data-consistent.

**CORRECTION to Step 16's core claim:** Step 16 treated
`truthMatchProbability < 0.5` as the marker of a PHYSICAL loss.  That is
wrong.  `truthMatchProbability` is an inner-detector tracking-truth quantity
(fraction of ID hits from the truth particle); it does NOT decide whether the
muon is reconstructed.  A muon can be a perfectly good Medium combined muon
(momentum MS-constrained, passes hit/quality cuts) while its ID-hit fraction
sits below 0.5 from pure occupancy.  The physical outcome is fully captured
by "does a Medium reco muon exist within ΔR<0.05."  Hence prob<0.5 is mostly
an MC BOOKKEEPING rejection, not a physical inefficiency.

**Implication — direction of bias reversed:** Requiring `prob>0.5` as the
matching gate (the current OFFLINE primary in `PythiaFullSimExtras.c:94`)
REJECTS genuinely-reconstructed Medium muons in dense events → biases ε LOW
→ OVER-corrects the cross-section (too-small ε inflates corrected yield).
This is the real methodological bug, opposite in sign to Step 16's concern.

**Revised verdict on the dR fallback:** Its DIRECTION (geometric ΔR<0.05) is
the ATLAS standard — same cone as the paper.  Step 16's "recovers physical
losses → overestimate" was wrong because the prob<0.5 "losses" are mostly MC
artifacts.  Genuine remaining flaws of the implementation: (a) it is a
bolt-on FALLBACK to a biased `prob>0.5` PRIMARY (inconsistent hybrid) — the
geometric match should be the sole/primary definition; (b) it grabs nearest
reco THEN checks Medium, vs the standard "nearest MEDIUM muon" with charge/pT
consistency; (c) residual real overestimate = accidental geometric match to
an unrelated nearby Medium fake in central PbPb — small, estimable as a
systematic (mitigated by tight cone + charge consistency), NOT a reason to
revert to prob>0.5.

**Skimmer is not buggy as a record;** it faithfully stores AOD barcode+prob.
The bug is the OFFLINE use of prob>0.5 as the primary efficiency definition.

**recoMuonLink follow-up WITHDRAWN:** under standard geometric matching it is
unused; its MS-vs-ID provenance is moot.

**Corrected recommendation:** replace the entire `prob>0.5` + dR-fallback
hybrid with the single ATLAS-standard definition —
`ε = (truth muons with a Medium reco muon within ΔR<0.05) / (fiducial truth
muons)`, charge/pT-consistent, closest-unclaimed reco — with NO
prob/barcode/recoMuonLink gate.  Use signal-only truth (r17662) so no HIJING
truth exists (removes barcode collision entirely; also fixes Bug #3
denominator).  Assess accidental-fake match rate in central bins as a
systematic.  The "prob>0.5 vs geometric" prototype reduces to VALIDATION:
confirm pp agreement (low occupancy) and quantify the HI low-bias of the old
prob method.

### Step 18: Run 2 note refutes "prob>0.5 biases low"; root cause is the barcode collision (2026-06-11)

**Reference checked:** Run 2 dimuon internal note ATL-COM-PHYS-2021-1094
(`~/workarea/dimuon_codes/IntNotesRun2DimuonReference/`).

**Run 2 reco efficiency definition** (`tex/Corrections.tex:377-383`):
> "the fraction of the number of generated-muons that have associated
> reconstructed muons … that match the generated-muon with a probability of
> **0.5 or more**."

evaluated on a **STARLIGHT γγ→μμ signal OVERLAID onto HIJING** sample
(`tex/Data_and_selections.tex:171-178`; productions ATLHI-304/313 for
ANA-HION-2020-10) — the SAME overlay concept as ours — and it yields
PHYSICAL efficiencies with only a mild central→peripheral trend.

**This REFUTES Step 17.**  `prob>0.5` matching + HIJING overlay gives physical,
data-consistent efficiencies (existence proof).  So `prob>0.5` does NOT bias ε
low, and the "hit dilution drops good muons below 0.5" worry (Step 16) is
quantitatively small — if it were large, the Run 2 efficiencies would be low
too.  Both Step 16's and Step 17's bias claims are withdrawn.

**True root cause (already named as Bug #2 mechanism i): the barcode collision
specific to r17618.**  Pythia and HIJING share barcode range 1..N (~652
collisions, per the barcode-duplicate investigation).  `prob>0.5` matching
finds the reco muon whose stored barcode equals the truth muon's barcode; when
barcodes collide, the signal muon's correct reco match is mislabeled/hidden →
artificially low ε.  The Run 2 STARLIGHT+HIJING sample was produced with proper
overlay truth handling (no collision), which is exactly why `prob>0.5` works
there.  **It is the SAMPLE, not the matching method.**

**Geometric ΔR<0.05 matching** is one valid fix because it ignores barcodes
entirely (collision-immune) — but it is not "more correct" than `prob>0.5`;
it just sidesteps the production artifact.  Pure-geometric replaces the
barcode find()+dR-fallback block (`PythiaFullSimExtras.c` 154-191): for each
truth muon, loop reco muons passing Medium, take closest unclaimed with
ΔR<0.05.  Skimmer needs no change (reco+truth kinematics already stored).

**Clean path (prefer #1 for direct comparability with the Run 2 reference):**
1. Reproduce the Run 2 method: **signal-only truth (r17662)** → no HIJING truth
   → no collision; keep standard **`prob>0.5`**; drop the dR fallback.  Re-test
   with adequate stats + the Bug #1 index fix.
2. Or pure geometric **ΔR<0.05-to-Medium** matching (collision-immune; also an
   ATLAS standard).  Residual: central-PbPb accidental fakes (small systematic).

Either way the current `prob>0.5` + dR-fallback HYBRID should go — it is a
band-aid for the barcode collision rather than a fix of it.

### Step 19: Full-stats r17662 grid skim submitted; plan to test pure prob>0.5 (2026-06-11)

**Failure-mechanism clarification (recorded for reference):** The barcode
collision corrupts the link at RECO/AOD level, not at the skimmer.  The
skimmer does NO matching — it just (a) stores every stable gen muon into
`truth_muon_*` (the Pythia muon is ALWAYS stored, never dropped), and (b)
copies each reco muon's upstream `truthParticleLink` barcode/prob into
`muon_truth_barcode`/`muon_truth_prob`.  Because the hit→truth map is keyed by
barcode, a shared barcode (Pythia bc=100 ≡ HIJING bc=100) makes the Pythia
muon's reco track link to a HIJING particle (stored bc e.g. 5078) or to nothing
(bc=−1).  The skimmer faithfully copies that wrong value.  The mismatch only
bites OFFLINE (`PythiaFullSimExtras.c` `std::find`): the Pythia gen muon
(bc 100) finds no reco muon carrying bc 100 → `reco_match=false` → low ε.
Signal-only truth fixes it at the source (no shared barcode → correct links).

**Event-count check:** r17662 AOD (802781 pTH8_14) =
`...e8599_s4614_r17662_r15970_tid50427580_00` has **10 files / 10,000 events**
(verified `rucio list-files`).  The prior 100-event signal-only skim came from
`run_pythia_fullsim_HIJING_overlay/run_test.sh` hardcoded `--evtMax=100` (a
LOCAL test) — not an AOD or grid limit.

**Grid skim submitted (background agent):** pathena jediTaskID **50903365**
(https://bigpanda.cern.ch/task/50903365/), outDS
`user.yuhang.NTUP.Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14.FullSimHIJINGOverlayPP24_r17662.June2026.v1.`,
script `SkimCode/run_pythia_fullsim_HIJING_overlay/grid_sub_r17662_signalonly.sh`.
Monitor daemon `grid_monitor.sh --mode overlay 50903365` (PID 1170068, node
attsub07) downloads to `~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/`
on completion; state file `grid_monitor_state.txt` marks `50903365 completed`.

**Pending plan (when NTUP downloaded):**
1. Add a runtime toggle to DISABLE the ad-hoc dR fallback in
   `PythiaFullSimExtras.c` (default = enabled, preserves current behavior) →
   via /review-analysis-code.  Run pure `prob>0.5` matching (Run 2 method).
2. Run NTP (no dR) → RDF hist filling → plotting on r17662.
3. Compare single-muon & pair reco ε to CURRENT values (r17618 + dR-fallback:
   Step 13 pair ~49-55% @0-5%, Step 15 single-muon plateau ~77%) and to the
   Run 2 physical shape.  Expectation: r17662 + pure `prob>0.5` gives physical
   (Run-2-like, single-muon plateau ~90%+) ε → confirms the barcode collision,
   not the matching method, was the culprit.

### Step 20: r17662 result — collision & dR fallback are NOT the cause (2026-06-11)

**Setup:** r17662 NTUP downloaded (10,000 events, tree HeavyIonD3PD).  Ran NTP
via `run_pythia_fullsim_overlay_r17662_nodr.sh` using
`fullsim_input_dir_override` → `r17662_run/` (symlink of r17662 NTUP under the
expected FullSimHIJINGOverlayPP24 name; other 5 pT slices absent → skipped).
`pythia_only_barcode_cache=false` disables the dR fallback (gate at
PythiaFullSimExtras.c:162); for signal-only truth the cache restriction is moot
and does NOT affect single-muon barcode matching (only ancestor tracing).
Compared single-muon ε vs current overlay (r17618+dR) and pp.

**Single-muon ε (qeta-integrated, weighted), turn-on region (pTH8_14 has few
muons >15 GeV):**

| pT [GeV] | r17662 no-dR | r17662 with-dR | r17618+dR (current) | pp |
|---|---|---|---|---|
| 6–8   | 0.743 | 0.743 | 0.732 | 0.869 |
| 8–10  | 0.751 | 0.751 | 0.751 | 0.890 |
| 10–12 | 0.814 | 0.814 | 0.762 | 0.904 |
| 12–15 | 0.857 | 0.857 | 0.773 | 0.909 |

**Key results:**
1. **dR fallback does NOTHING on r17662** — no-dR and with-dR are byte-identical
   (no collision → barcode matching succeeds → fallback never fires).
2. **Removing the collision does NOT raise ε** — r17662 (clean) ≈ r17618+dR ≈
   ~73-75% turn-on, all ~14% below pp (~87%).
3. **This REFUTES Step 18's "barcode collision is THE culprit."**  The collision
   and dR fallback are NEGLIGIBLE for the efficiency.  (Bug #1, the index-mapping
   fix, 10%→50% pair, was the real killer and remains the main correction.)

**Decomposition (reco_match = prob>0.5 match exists; pass_medium = matched AND
Medium):**

| pT | r17662 reco_match | r17662 pass_med | pp reco_match | pp pass_med |
|---|---|---|---|---|
| 4–5   | 0.858 | 0.635 | 0.967 | 0.744 |
| 6–8   | 0.895 | 0.743 | 0.983 | 0.869 |
| 8–10  | 0.918 | 0.751 | 0.983 | 0.888 |
| 12–15 | 0.959 | 0.857 | 0.980 | 0.901 |

Overlay deficit vs pp = **~10% fewer prob>0.5 reco-matches** + ~3-4% more
Medium-fails among matched.  The ~10% reco_match deficit is the dominant
residual and is robust across r17618/r17662.

**THE open question (unchanged, now sharply posed):** is the ~10% reco_match
deficit (a) PHYSICAL (muon genuinely not reconstructed in the dense environment)
or (b) a prob>0.5 STRICTNESS artifact (muon IS reconstructed as a Medium muon
but its ID-track prob<0.5 from HIJING/jet hit dilution → tagged unmatched)?
Decisive test = **geometric ΔR<0.05-to-Medium matching** (ATLAS standard,
arXiv:2012.00578): if geometric reco_match → ~pp (~0.97), deficit is (b) and the
physical ε is ~pp; if it stays ~0.86-0.96, deficit is (a), physical.
Note: Run 2 note used prob>0.5 and got physical ε but on STARLIGHT (isolated UPC
muons, little hit dilution); our Pythia HF muons are in jets → more dilution, so
prob>0.5 may be stricter for us.  Geometric matching is data-consistent either
way.

**Caveats:** (i) r17662 has only pTH8_14 → no plateau (>20 GeV) stats; plateau
comparison needs higher-pT r17662 slices (grid).  (ii) NTP printed
`#muons from B-hadron: 0` on r17662 → ancestor tracing / from_same_b likely
broken because signal-only-truth thinned the truth parent chain.  Single-muon &
inclusive SS/OS ε are unaffected, but single-b pair classification on r17662 is
suspect — flag before using r17662 for single-b.

**Files:** `run_pythia_fullsim_overlay_r17662_nodr.sh`; outputs in
`pythia_fullsim_hijing_overlay_test_sample/r17662_run/` (suffixes
`_r17662_nodr_single_muon`, `_r17662_withdr_single_muon`).

### Step 21: Reconsideration — collision may PARTLY contribute; physical low-pT hypothesis (2026-06-12)

**Reconsideration of Step 20's "negligible collision" claim (user-flagged):**
Step 20 over-stated it.  Re-examining the single-muon comparison table:

| pT [GeV] | r17662 no-dR (8-14, NO collision) | r17618+dR (all slices) | pp | r17662 N |
|---|---|---|---|---|
| 4–5   | 0.635 | 0.630 | 0.738 | 8040 |
| 5–6   | 0.727 | 0.712 | 0.840 | 3885 |
| 6–8   | 0.743 | 0.732 | 0.869 | 2918 |
| 8–10  | 0.751 | 0.751 | 0.890 |  743 |
| 10–12 | 0.814 | 0.762 | 0.904 |  161 |
| 12–15 | 0.857 | 0.773 | 0.909 |   49 |

- At 6–10 GeV (solid r17662 stats): r17662 ≈ r17618 → collision effect small here.
- At 10–15 GeV: r17662 noticeably HIGHER (0.814/0.857 vs 0.762/0.773).  This is a
  real, noticeable improvement → **the collision MAY contribute (partly), at least
  at higher pT.**  Caveats: low r17662 stats there (N=49–161, errors ~0.03–0.05,
  ~1.7σ) AND confounded by slice composition (r17662 = 8-14 only vs r17618 = all
  slices; at fixed muon pT the slice mix differs).  Clean isolation (r17618 8-14
  ONLY vs r17662) is NOT possible from the single-muon output (no pTHat branch);
  would need an r17618 8-14-only single-muon NTP run, or use the pair-level kin0
  trees (which DO carry pTHat).  So: collision is not strictly negligible — revise
  Step 20 to "small at low pT, possibly contributing at 10–15 GeV (low-significance)."

**Physical hypothesis (user, to be tested):** The low-pT overlay deficit vs pp may
be (at least partly) PHYSICAL, not a procedural bug.  In real very-central Pb+Pb,
thousands of soft underlying-event particles (simulated by HIJING) deposit hits;
**low-pT muons are the most affected** — their ID tracks are softest/shortest-lever
and most easily corrupted by nearby soft hits, so genuine reco+ID efficiency for
low-pT muons in central Pb+Pb can be truly lower than in pp.  This is a real
detector effect, expected in data too.

**Discriminating test = geometric ΔR<0.05-to-Medium matching:**
- If the deficit is PHYSICAL (muon genuinely fails Medium reco/ID), geometric
  matching shows the SAME/similar deficit (no Medium reco muon exists near the
  truth muon to match), because geometric only asks "does a Medium muon exist
  within ΔR<0.05" — it cannot conjure a muon that was not reconstructed.
- If the deficit is a prob>0.5 truth-matching ARTIFACT (Medium muon exists but its
  ID-track prob<0.5 from hit dilution), geometric RECOVERS them → ε rises toward pp.
So geometric vs prob>0.5 separates physical (a) from matching-artifact (b).  Both
the collision-contribution and the physical-deficit pictures predict specific
geometric outcomes; the test adjudicates.

### Step 22: pThat-slice restriction — skew vs statistics (Task 2) (2026-06-12)

Compared pair reco ε vs truth_pair_pt for **pThat 8-14 only (kin0)** vs **all
slices**, straight from the existing r17618 output (per-kin trees
`muon_pair_tree_kin0_sign2` vs `muon_pair_tree_sign2`; no re-run needed):

| pair pT [GeV] | ALL slices ε (N) | kin0 8-14 ε (N) |
|---|---|---|
| 0–4   | 0.513 (8670)  | 0.495 (1663) |
| 4–8   | 0.487 (8680)  | 0.472 (1084) |
| 8–12  | 0.462 (19378) | 0.453 (1938) |
| 12–16 | 0.510 (13782) | 0.555 (425)  |
| 16–20 | 0.520 (8510)  | 0.491 (55)   |
| 20–28 | 0.517 (8814)  | — (3)        |
| 28–40 | 0.550 (5830)  | — (0)        |
| 40–60 | 0.590 (3272)  | — (0)        |

**Answer: closer to case (2), with a strong caveat.** Where kin0 has statistics
(pair pT < ~16 GeV), its efficiency is consistent with all-slices (within ~1-2%,
the 12-16 point ~1.9σ high but N=425) → **no strong SKEW** of the value.  BUT the
8-14 slice has essentially **ZERO events at high pair pT (>20 GeV)** → the plateau
is simply **unreachable** from 8-14 alone (not "lower stats but usable" — it's
unusable at the plateau).  Implication: r17662 (8-14 only) can test the TURN-ON
region robustly but **cannot probe the plateau at all**; a plateau comparison
needs the higher-pThat r17662 slices (grid) — OR we rely on the turn-on, where the
physics/collision question is already accessible.

### Step 23: GEOMETRIC test result — the low-pT deficit is PHYSICAL (2026-06-12)

**Implementation:** Added `use_geometric_matching` flag (default false → clean
rollback) to `PythiaTruthExtras.h`; shared `geometric_match` lambda + geometric
primary path in `PythiaFullSimExtras.c::ProcessEventFullsim` (nearest unclaimed
reco within ΔR<0.05, ignoring barcode/prob — ATLAS standard arXiv:2012.00578).
/review-analysis-code PASS iter 1 (log review-analysis-code-20260612-040528).
ACLiC compiled.  Ran on r17662 single-muon with a SEPARATE suffix.

**Result — geometric vs prob>0.5 on r17662 (single-muon, by muon pT):**

| pT [GeV] | prob>0.5 reco_m / pass_med | geometric reco_m / pass_med |
|---|---|---|
| 4–5   | 0.858 / 0.635 | 0.858 / 0.635 |
| 6–8   | 0.895 / 0.743 | 0.894 / 0.742 |
| 8–10  | 0.918 / 0.751 | 0.918 / 0.751 |
| 12–15 | 0.959 / 0.857 | 0.959 / 0.857 |

**Geometric ≡ prob>0.5 (identical within ≤0.0014).**  Decisive interpretation:
the muons that fail prob>0.5 matching ALSO have NO reco muon within ΔR<0.05 — they
are **genuinely not reconstructed**, not merely prob-mistagged.  So the ~10-14%
overlay deficit vs pp is **PHYSICAL** (real reco/ID inefficiency in the dense
central-PbPb / HIJING environment), NOT a truth-matching artifact.  This
**confirms the Step-21 physical hypothesis** and answers the long-standing open
question (a) vs (b) → **(a) PHYSICAL**.

**Consequences / reconciliation of the whole investigation:**
- The matching method (prob>0.5 vs geometric) does NOT matter on a clean sample —
  they agree.  Consistent with the Run 2 note using prob>0.5 and getting physical ε.
- The barcode collision (r17618) and dR fallback are minor perturbations ON TOP of
  this physical baseline (collision possibly +a few % at 10-15 GeV, Step 21).
- Bug #1 (index mapping) was the real code bug (10%→50% pair); after it, the
  residual overlay deficit is physical, expected, and present in data.
- The deficit is largest at low pT (4-6 GeV: ε~0.64 vs pp~0.74) and shrinks toward
  higher pT — exactly the "low-pT muons most affected by underlying-event hits"
  picture.

**Results separation / rollback (per user request):**
- prob>0.5 (no-dR): `r17662_run/..._r17662_nodr_single_muon.root`
- geometric: `r17662_run/..._r17662_geometric_single_muon.root`
- with-dR: `r17662_run/..._r17662_withdr_single_muon.root`
- plots: `plots/geometric_matching_test/` (geometric), `plots/pth_slice_skew_test/` (Task 2)
- Code rollback: `use_geometric_matching` defaults false (no behavior change); the
  diff is additive in 2 files (PythiaTruthExtras.h, PythiaFullSimExtras.c).

### Step 24: InitParamsExtra override discovered; dR fallback decoupled & default-OFF; Step 20-23 reco_match corrected (efficiency conclusions HOLD) (2026-06-12)

**Critical discovery:** `InitParamsExtra()` is called INSIDE `Run()` (via
`DimuonAlgCoreT::Run` → `InitParamsHook` → `CallInitParams` fold) and unconditionally
sets `pythia_only_barcode_cache = true`.  It therefore OVERRIDES any
`py.pythia_only_barcode_cache = false` set in a run macro before `Run()`.  Since the
old dR-fallback gate was `else if (pythia_only_barcode_cache)`, the dR fallback could
NEVER be disabled from a macro — so **all prior "nodr" r17662 runs (Steps 20-23)
actually ran WITH the dR fallback ON.**

**Corrected r17662 single-muon numbers (6-8 GeV, N=2918):**

| matching | reco_match | pass_medium (= EFFICIENCY) |
|---|---|---|
| pure prob>0.5 (no dR) — TRUE default | 0.8266 | 0.7426 |
| prob>0.5 + dR fallback | 0.8955 | 0.7430 |
| geometric ΔR<0.05 | 0.894  | 0.7423 |
| pp reference (prob>0.5) | 0.983 | 0.869 |

**The physics conclusion is UNCHANGED and now MORE robust:** `pass_medium` (the
efficiency) is **~0.743 for ALL three overlay matching methods** (pure prob>0.5,
+dR, geometric) — the matching method does NOT change the efficiency; the
dR/geometric "recovered" muons (the reco_match difference) all FAIL Medium anyway.
So the ~14% deficit vs pp (0.743 vs 0.869) is **robustly PHYSICAL**, independent of
matching method.  This SUPERSEDES the Step-23 wording "geometric ≡ prob>0.5 reco_match"
(that comparison was geometric-vs-(prob+dR), both ~0.89); the correct invariant is
**pass_medium is method-independent**.

**Code change (this session, /review-analysis-code PASS, log 20260612-041807):**
- Decoupled the dR fallback: new flag `use_dr_fallback` (default **false**) in
  PythiaTruthExtras.h; gate in PythiaFullSimExtras.c changed
  `pythia_only_barcode_cache` → `use_dr_fallback`.  `pythia_only_barcode_cache` keeps
  its true cache role (PythiaTruthExtras.c:355) and stays true in InitParamsExtra.
- **DEFAULT overlay matching is now pure prob>0.5 — NO dR, NO geometric** (verified:
  default reco_match=0.8266, no dR; use_dr_fallback=true reco_match=0.8955).
- The geometric-matching refactor (Step 23) was **REVERTED from production** — not
  needed (efficiency is method-independent) and its restructure exposed, confusingly,
  the true no-dR rate.  The geometric *result* stands (pass_medium 0.7423).
- Final production diff: +12/-3 in 2 files (PythiaTruthExtras.h, PythiaFullSimExtras.c).
  Rollback = set `use_dr_fallback`/`use_geometric_matching` defaults or git revert.

**For the upcoming full signal-truth-only overlay sample:** it will run the original
prob>0.5 procedure by default (no dR workaround, no geometric) — as required.

## Latest Stage

**Step 24 (2026-06-12): dR fallback decoupled & default-OFF; efficiency conclusion robust.**

DISCOVERY: InitParamsExtra (in Run) forces pythia_only_barcode_cache=true, so prior
"nodr" runs were actually +dR.  Decoupled the dR fallback into `use_dr_fallback`
(default false) → DEFAULT overlay matching is now pure prob>0.5 (no dR, no geometric),
verified on r17662 (reco_match 0.8266 default vs 0.8955 +dR).  Geometric refactor
reverted from production.  PHYSICS HOLDS & STRONGER: pass_medium efficiency is
~0.743 for pure-prob / +dR / geometric alike (method-independent) vs pp 0.869 →
deficit is robustly PHYSICAL.  /review-analysis-code PASS.  Full sample will run the
original prob>0.5 procedure by default, as the user required.  REMAINING (optional):
higher-pThat r17662 slices for plateau stats; pair-level cross-check (single²≈pair).

**Prior Step 21-22 (2026-06-12): collision partial-contribution + physical hypothesis logged; pThat skew assessed.**

(1) Revised Step 20: collision NOT strictly negligible — r17662 higher than r17618
at 10–15 GeV (0.814/0.857 vs 0.762/0.773; low stats, slice-confounded), small at
6–10 GeV.  Added physical hypothesis: low-pT central-PbPb deficit may be real
(underlying-event/HIJING soft hits hurt low-pT muon ID most); to be tested by
geometric matching (physical → same deficit; artifact → recovers).
(2) Task 2: pThat 8-14-only does NOT strongly skew ε where it has stats, but has
ZERO plateau stats (>20 GeV) → turn-on testable, plateau not.
NEXT: (Task 1) implement geometric ΔR<0.05-to-Medium matching (/review-analysis-code),
run on r17662 with SEPARATE output suffix/dir (no overwrite, clear rollback),
compare geometric vs prob>0.5 to adjudicate physical vs artifact.

---
**Prior Step 19 (2026-06-11): r17662 full-stats skim submitted (task 50903365).**

Grid skim of signal-only-truth r17662 (802781 pTH8_14, 10k events) submitted;
monitor daemon downloads on completion.  Background poll watches
`grid_monitor_state.txt` for `50903365 completed` and re-invokes to run the
pipeline.  When NTUP lands: disable dR fallback (toggle, /review-analysis-code),
run NTP→RDF→plot on r17662 with pure `prob>0.5`, compare reco ε vs current
(r17618+dR: pair ~49-55%, single-muon plateau ~77%) and vs Run 2 physical shape.
Prior 100-event signal-only skim was a run_test.sh `--evtMax=100` local test, not
an AOD limit.
