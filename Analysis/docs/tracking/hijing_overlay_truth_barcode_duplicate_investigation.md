# HIJING Overlay: Truth / Barcode Duplicate Investigation

## Objective

Investigate truth-record and barcode-duplicate behaviour across the three
new overlay r-tags (r17618, r17662, r17663) to understand how
`StandardSignalOnlyTruth` and no-overlay configurations affect truth
content, and whether barcode collisions or truth contamination exist.

## Context

Pythia 5.36 TeV pp hQCD DiMu FullSim samples were reconstructed with
HIJING overlay under three different reco configurations.  Prior work
(see `hijing_overlay_debug_summaries.md`) established that the original
overlay AODs merge Pythia + HIJING truth into a single `TruthParticles`
container, with HIJING barcodes at bc > 200 000.  The new r-tags provide
two alternative truth strategies (signal-only truth, no overlay at all)
for comparison.

## Datasets & Skimming

### R-tags

All three use **Athena 24.0.58**, conditions `OFLCOND-MC23-SDR-RUN3-05`,
geometry `ATLAS-R3S-2021-03-02-00`, HI trigger menu
`Dev_HI_run3_v1_TriggerValidation_prescale`, `ConditionsRunNumber=460000`,
beamspot `posZ=-3.3 mm`.  Created by Ewelina for ATLHI-576 (impact
parameter 0–5).

Full AMI dumps saved at:
`/gpfs/mnt/atlasgpfs01/usatlas/data/yuhanguo/pythia_fullsim_hijing_overlay_test_sample/r17618.txt`,
`r17662.txt`, `r17663.txt`.

| R-tag | Description | Overlay | digiSteeringConf | inputCavernHitsFile |
|-------|-------------|---------|------------------|---------------------|
| r17618 | Clone of r17038, beamspot fix | PileUp=True | — | Hijing UCC e8613_s4614_s4258 |
| r17662 | Clone of r17618 + StandardSignalOnlyTruth | PileUp=True | StandardSignalOnlyTruth | Hijing UCC e8613_s4614_s4258 |
| r17663 | Clone of r17618 w/o HIJING overlay | PileUp=False | — | (none) |

Key difference between r17618 and r17662: `StandardSignalOnlyTruth`
strips HIJING truth particles from the output truth containers, keeping
only signal-generator (Pythia) truth records.  Detector-level overlay
effects are still present (PileUp=True), so reco tracks, calorimeter
deposits, and HIEventShape all include HIJING contributions.

Key difference for r17663: no overlay at all — `PileUp=False`,
`DoXingByXingPileUp=False`, no `inputCavernHitsFile`.  This is pure
signal-only reconstruction.

### AOD datasets (pTH8–14 slice)

| R-tag | Dataset | Files | Size/file | Replica |
|-------|---------|-------|-----------|---------|
| r17618 | `mc23_5p36TeV.802781.…e8599_s4614_r17618_r15970_tid50035726_00` | 10 | ~8 GB | TOKYO-LCG2_DATADISK |
| r17662 | `mc23_5p36TeV.802781.…e8599_s4614_r17662_r15970_tid50427580_00` | 10 | ~6 GB | AGLT2_DATADISK |
| r17663 | `mc23_5p36TeV.802781.…e8599_s4614_r17663_r15970_tid50446417_00` | 10 | ~500 MB | INFN-T1_DATADISK |

Note the dramatic size difference: r17663 (no overlay) files are ~500 MB
vs ~6–8 GB for overlay samples.

### Skimming procedure

Skimmer: `SkimCode/run_pythia_fullsim_HIJING_overlay/TrigRates_CA.py`
with `TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24`.

Build: AthAnalysis 25.2.89, built via `SkimCode/setup_25.sh`
(`acmSetup` + `acm compile` in `build_25/`).

Run mode settings (`do_pp_MC_fullsim_hioverlay24=True`):
- `is_HION=True` → TRT cut-off enabled on muon selection
- `is_MC=True` → trigger processing off, truth storage on
- `RunYear=2024` → 2024 HIEventShape centrality calibration
- `StoreZdc=0` → ZDC disabled (not present in overlay MC AODs)
- `StoreTruth=15` (bits 1+2+4+8) → full truth with parents & links
- `HIEventShapeContainerKey="HIEventShape"` → reads HI event shape

#### Test skimming results (100 events each)

Run via `run_test.sh` helper script.  All three succeeded (exit code 0).

| R-tag | Suffix | Output | Output size |
|-------|--------|--------|-------------|
| r17618 | `_orig` | `Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14._orig.root` | 213 MB |
| r17662 | `_signalOnlyTruth` | `Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14._signalOnlyTruth.root` | 1.6 MB |
| r17663 | `_no_overlay` | `Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14._no_overlay.root` | 1.7 MB |

Output location:
`/gpfs/mnt/atlasgpfs01/usatlas/data/yuhanguo/pythia_fullsim_hijing_overlay_test_sample/`

#### Skimming gotchas

1. **Athena environment in non-interactive shells:** `acmSetup` and
   `asetup` do not fully propagate the Athena runtime PATH in
   non-interactive bash.  The AthAnalysis release `bin/` directory
   (containing `athena.py`) is not added to PATH by the local WorkDir
   `setup.sh`.  Workaround: source the saved setup
   (`build_25/.asetup.save`) then the local `setup.sh`, or explicitly
   prepend the release bin:
   `export PATH=/cvmfs/.../AthAnalysis/25.2.89/InstallArea/x86_64-el9-gcc14-opt/bin:$PATH`.

2. **`libexcabort.so` not found:** The `athena.py` wrapper script
   (`USEEXCABORT=1` by default) looks for `libexcabort.so` in
   `LD_LIBRARY_PATH`.  Not present in AthAnalysis 25.2.89 non-interactive
   setup.  Fix: pass `--no-excabort` to `athena.py`.

3. **Remote xrootd I/O:** Overlay AODs are 6–8 GB each and only
   available at remote sites (Tokyo, AGLT2, INFN).  INFN xrootd
   (`xrootd-atlas.cr.cnaf.infn.it`) was completely unresponsive for
   metadata peeking.  Solution: `rucio download` the file locally first
   (only ~500 MB for r17663, downloads in ~13 s).

4. **HIEventShape survives without HIJING overlay:** r17663 (no overlay,
   `PileUp=False`) still produces `HIEventShape` because the reco tag
   includes `preInclude: HIRecConfig.HIModeFlags.HImode`.  The skimmer's
   `EventShape::Process()` retrieves the container without error.  FCal_Et
   will be tiny (signal-only calorimeter deposits), and centrality will be
   meaningless — but the code does **not** crash.  No need to set
   `HIEventShapeContainerKey=""` for r17663.

5. **Output size difference:** r17618 (_orig) produces 213 MB / 100 events
   because `StoreTruth=15` stores all ~70 k truth particles per event
   (Pythia + HIJING merged in `TruthParticles`) plus thousands of reco
   tracks from the HIJING underlying event.  r17662 and r17663 produce
   only ~1.6 MB because truth and reco track counts are dramatically
   smaller (signal-only).

6. **`meTrk Not Found for Combined muon`:** Harmless ERROR-level message
   from `TrigRates.cxx` when a combined muon has no ME track pointer.
   Present in all three samples.  Not a failure.

7. **MuonEfficiencyTool warnings:** r17663 emits warnings about
   `run number 999999` not matching any SF period.  This is because
   `ConditionsRunNumber=460000` sets the run to a dummy value, and the
   efficiency tool cannot find scale factors for it.  Harmless for
   truth-level studies.

## Sub-steps

1. Document datasets and skimming procedure ← **done**
2. Compare truth containers across the three NTUPs (TruthParticles count,
   barcode ranges, parent chains)
3. Check for barcode duplicates or collisions in r17618 (orig overlay)
4. Verify signal-only truth cleanliness in r17662 vs r17663
5. Quantify reco-level differences (track multiplicities, FCal_Et, centrality)
6. Summarize findings and recommend which r-tag to use for production

## Accumulated Findings

### Step 1: Dataset & skimming (2026-06-04)

- All three r-tags skim successfully with `TrigRates_CA.py` in
  `ppmcfullsim_hioverlay24` mode, no code changes needed.
- r17662 (StandardSignalOnlyTruth) dramatically reduces output size by
  stripping HIJING truth (1.6 MB vs 213 MB per 100 events).
- r17663 (no overlay) does NOT crash on HIEventShape — the container
  exists even without HIJING, produced by the HI reco preInclude.
- Full-AOD re-skim: 1000 events each. r17618=2.1 GB, r17662=17 MB,
  r17663=16 MB. No common eventNumbers across any pair of files (different
  EVNT inputs).

### Step 2: Truth tracing & category comparison (2026-06-04)

**Code change:** Added `fullsim_input_dir_override` to `PythiaAlgCoreT.h/.c`
(2-line change) so the analysis can read NTUPs from a custom directory
without modifying `FullSimSampleType.h`.

**Setup:** Created 3 temp subdirectories under
`pythia_fullsim_hijing_overlay_test_sample/` (`rtag_orig/`,
`rtag_signalOnlyTruth/`, `rtag_no_overlay/`), each with a symlink
`Pythia_5p36TeV_pp_hQCD_DiMu_pTH8_14.FullSimHIJINGOverlayPP24.NTUP.root`
→ the corresponding test NTUP.  Used `extra_output_suffix` to differentiate
outputs.

**NTP processing results (pTH8_14 only, 1000 events each):**

| R-tag | Suffix | B-hadron muons | Truncated | CPU time |
|-------|--------|----------------|-----------|----------|
| r17618 | `_rtag_orig` | 915 | 27 (2.95%) | 30 s |
| r17662 | `_rtag_signalOnlyTruth` | 832 | 0 | 0.78 s |
| r17663 | `_rtag_no_overlay` | 796 | 0 | 0.81 s |

Key: r17618 has 2.95% truncated B-hadron histories from HIJING barcodes
(bc 7154–61923) — these are HIJING-origin muons being truth-matched to
signal reco muons (or HIJING muons entering the truth-pair loop).  r17662
and r17663 have zero truncations.

**Pair counts (SS / OS):**
- r17618: 145 / 551
- r17662: 111 / 508
- r17663: 94 / 507

r17618 has more pairs due to HIJING-origin muons forming additional
truth-matched reco pairs.

**Category comparison plots** (8 plots with error bars in
`plots/truth_origin_comparison/rtag_*.png`, plotter:
`plot_rtag_comparison.C`):

**All three r-tags agree within statistical errors on all flavor, parent
group, and most origin categories.**  The only statistically significant
difference is the **"others" origin category**:

| Category | r17618 (orig) | r17662 (sigTruth) | r17663 (noOverlay) |
|----------|---------------|--------------------|--------------------|
| others (OS) | 0.058 ± 0.010 | 0.020 ± 0.006 | 0.014 ± 0.005 |
| others (SS) | 0.234 ± 0.040 | 0.063 ± 0.024 | 0.074 ± 0.028 |

- OS "others": r17618 is ~3σ above r17662/r17663 (which agree with each
  other).
- SS "others": r17618 is ~3.7σ above r17662/r17663.
- The "others" excess in r17618 comes at the expense of "fc" (flavor
  creation): SS fc = 0.276 ± 0.044 (orig) vs 0.450 ± 0.064 (sigTruth),
  ~2.2σ pull.  Pairs whose ancestor tracing fails (due to HIJING truth
  contamination) get dumped into "others" instead of being classified as fc.
- r17662 and r17663 agree with each other in all bins.
- All other categories (flavor, parent group) agree within 1σ across
  all three r-tags.

**Physical interpretation:** The 2.95% truncated B-hadron histories in
r17618 (HIJING barcodes bc 7154–61923) cause ancestor tracing failures
that inflate the "others" origin category, particularly for SS pairs.
r17662 and r17663 eliminate this contamination entirely.  Since r17662
preserves detector-level HIJING overlay effects in reco while r17663
does not, r17662 is the preferred choice for production.

## Temporary Code Changes (to be undone after testing)

*(Track all temporary modifications here. Remove each entry once undone.
If we adopt r17662 for production, none of these changes need to persist.)*

| File | Change | Purpose | Status |
|------|--------|---------|--------|
| `PythiaTruthExtras.h` | Added `bool pythia_only_barcode_cache = false;` member | Step 3a: flag to restrict barcode cache | **done** |
| `PythiaTruthExtras.c:GetParticleIndex` | When `pythia_only_barcode_cache`, cache only up to first bc > 200k | Step 3a: test barcode duplication hypothesis | **done** |
| `run_rtag_comparison.sh` | Added optional 3rd arg `pythia_only_bc` | Step 3a: pass flag to NTP processing | **done** |

## Ruled Out

(none yet)

## Latest Stage

**Step 3 complete. (2026-06-04)**

### Step 3a: Barcode duplication in r17618

**Truth container structure in r17618 (event 0, 81030 particles):**

```
Pythia truth:    idx 0..650    (bc 1..651)
Pythia Geant4:   idx 651..662  (bc 200001..200171)
HIJING truth:    idx 663..76522 (bc 1..75860, RESTARTS FROM 1!)
HIJING Geant4:   idx 76523..81029 (bc 200027+)
```

**652 duplicate barcodes** (barcodes 1–651 all appear twice: Pythia idx +
HIJING idx).  `barcode_to_index_cache.emplace()` stores the **first**
occurrence (Pythia), so Pythia muons trace correctly.  But:

**Reco muon truth-matching reveals the mechanism:**
- Reco muons truth-matched to Pythia particles (e.g. bc=633) have
  ambiguous barcodes but resolve to the Pythia index (first occurrence) → OK
- Reco muons truth-matched to HIJING particles (bc=67878, 11699, etc.
  at HIJING+ indices) are unique → they enter truth-pair analysis but
  their ancestor chains are HIJING particles, causing truncation.
- ~4.7 reco muons/event in r17618 (due to HIJING overlay), many from HIJING.

**Pythia-only barcode cache fix (temporary):**
Added `pythia_only_barcode_cache` flag to `PythiaTruthExtras`.  When
set, `GetParticleIndex` caches only Pythia truth (indices before first
bc > 200k).  HIJING muons now fail to resolve in `SingleMuonAncestorTracing`
(SAT warning + skip), instead of silently misclassifying.

**Results:** Reran r17618 with `pythia_only_barcode_cache=1`:
- **0 truncated B-hadron histories** (was 27 = 2.95%)
- 888 B-hadron muons (was 915; 27 fewer = HIJING muons no longer entering B-hadron tracing)
- Many SAT warnings for HIJING muon barcodes not in cache → correctly skipped.

**Before/After comparison** (8 plots at `before_after_*.png`):

| Category | Before (all bc) | After (Pyth-only) | r17662 | r17663 |
|----------|-----------------|---------------------|--------|--------|
| others (OS) | 0.058 ± 0.010 | **0.015 ± 0.005** | 0.020 ± 0.006 | 0.014 ± 0.005 |
| others (SS) | 0.234 ± 0.040 | **0.103 ± 0.027** | 0.063 ± 0.024 | 0.074 ± 0.028 |
| not_both_open_HF (OS) | 0.479 ± 0.029 | 0.523 ± 0.031 | 0.461 ± 0.030 | 0.462 ± 0.030 |
| not_both_open_HF (SS) | 0.117 ± 0.028 | 0.248 ± 0.041 | 0.081 ± 0.027 | 0.096 ± 0.032 |

- OS "others" now agrees with r17662/r17663 (0.015 vs 0.020/0.014).
- SS "others" improved from 0.234 to 0.103 (now ~1σ from r17662/r17663).
- `not_both_open_HF` increases because HIJING muons in pairs (whose
  barcode can't be resolved) now get `parent_group = -10` → they fall
  into `not_both_open_HF` instead of being misclassified.

**Conclusion (3a):** Barcode duplication is the primary root cause of the
"others" excess in r17618.  The fix (Pythia-only barcode cache)
eliminates the misclassification for OS and dramatically reduces it for SS.
The remaining SS discrepancy comes from the higher pair count in r17618
(more HIJING-origin reco muons forming pairs).

### Step 3b: r17662 vs r17663 disagreement

**Investigation findings:**

1. **Different EVNT inputs confirmed:** r17662 eventNumbers start at
   2601,2646,2640,... while r17663 starts at 1,7,11,...  These are
   different MC draws → statistical fluctuations are expected.

2. **No barcode duplication in either:** r17662 has 0 duplicates (331
   truth particles/event, bc 1..319 + Geant4 bc > 200k).  r17663 also
   0 duplicates (486 particles/event, bc 1..475 + Geant4 bc > 200k).
   r17663 has more truth particles because the full Pythia shower is
   preserved (r17662's StandardSignalOnlyTruth may trim some particles).

3. **Reco muon multiplicity differs drastically:**
   - r17662: avg **4.7 reco muons/event** (HIJING detector overlay →
     many HIJING tracks reconstructed as muon candidates)
   - r17663: avg **1.9 reco muons/event** (pure signal → only Pythia muons)
   - This is the primary driver of the pair count difference
     (r17662=111SS/508OS vs r17663=94SS/507OS).

4. **Truth content is equivalent:** Both have Pythia-only truth with no
   HIJING contamination.  The truth-level ancestor tracing works
   identically for both (0 truncated chains each).

**Conclusion (3b):** r17662 and r17663 disagreements are purely
statistical (different EVNT inputs, 1000 events) and driven by the
reco muon selection difference (HIJING detector overlay present in r17662
but not r17663).  There is **no bug** in the truth analysis for these
two samples.  All categories agree within 1–2σ, consistent with
statistical expectations for ~100–500 pairs.

### Summary of findings

The truth analysis code has **no bug for signal-only truth samples**
(r17662, r17663).  The disagreement between r17618 (original) and the
other two is entirely caused by **HIJING truth barcode duplication**:
Pythia and HIJING truth particles share barcodes 1..N, and HIJING muons
entering the truth-pair analysis get misclassified because their ancestor
chains lead into HIJING territory that the tracing code was not designed
to handle.

**Suggested fix for production:**
- Use **r17662 (StandardSignalOnlyTruth)** which eliminates HIJING truth
  entirely while preserving detector-level overlay effects.
- No code change needed — the existing truth analysis works correctly
  with r17662 out of the box.
- The `pythia_only_barcode_cache` fix is a workaround for r17618 and
  should be removed once r17662 is adopted.

**Files:**
- Backups: `rtag_orig/*.bak_before_pythia_only`, `plots/bak_before_pythia_only/`
- New output: `rtag_orig/*_rtag_orig_pythiaonly.root`
- Plots: `plots/truth_origin_comparison/before_after_*.png` (8 plots)
