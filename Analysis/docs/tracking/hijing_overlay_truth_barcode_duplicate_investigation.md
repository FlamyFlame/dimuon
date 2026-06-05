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

## Ruled Out

(none yet)

## Latest Stage

**Step 2 complete.**  All three r-tags processed, 8 comparison plots
produced.  Next: deeper truth container comparison (sub-step 2–4) or
summarize findings and recommend r-tag for production.
