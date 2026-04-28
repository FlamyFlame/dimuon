# Pythia Fullsim Overlay Analysis

Pythia di-muon pairs from open heavy-flavour decays, fully simulated and reconstructed with HIJING (or Zmumu / data) overlay. Measures muon-pair reconstruction efficiency and detector response as a function of centrality in Pb+Pb conditions.

## Sample types

Controlled by `FullSimSampleType` enum (`MuonObjectsParamsAndHelpers/FullSimSampleType.h`):

| Enum value | Label | Overlay | Input directory (test) |
|-----------|-------|---------|----------------------|
| `hijing` | `hijing_overlay_pp24` | HIJING 0–5 fm impact parameter (ultra-central) | `usatlasdata/pythia_fullsim_hijing_overlay_test_sample/` |
| `zmumu` | `zmumu_overlay_pp24` | Z->mumu | `usatlasdata/pythia_fullsim_zmumu_overlay_test_sample/` |
| `data` | `data_overlay_pp24` | Real data | `usatlasdata/pythia_fullsim_data_overlay_test_sample/` |

Helper functions `FullSimSampleInputDir()`, `FullSimSampleLabel()`, `FullSimSamplePlotDir()`, etc. derive all I/O paths from the enum.

**Current limitation:** All paths point to test samples (6 NTUP files, ~10k events each). When full overlay samples become available, update `FullSimSampleType.h` paths and the Condor `.sub` file queue count.

## Key files

### NTuple processing (Stage 1)

| File | Role |
|------|------|
| `NTupleProcessingCode/PythiaAnalysisClasses.h` | Defines `PythiaFullSimOverlayAnalysis` (concrete class) |
| `NTupleProcessingCode/PythiaAlgCoreT.{h,c}` | Template algorithm core; drives event loop |
| `NTupleProcessingCode/PythiaFullSimExtras.{h,c}` | Reco-truth matching, muon quality cuts |
| `NTupleProcessingCode/PythiaFullSimOverlayExtras.{h,c}` | Overlay-specific: centrality, FCal_Et, beam type, weight normalization |
| `NTupleProcessingCode/PythiaTruthExtras.{h,c}` | Truth ancestry tracing (HF origin tagging) |
| `NTupleProcessingCode/run_pythia_fullsim_overlay_condor.sh` | Condor worker script |
| `NTupleProcessingCode/run_pythia_fullsim_overlay.sub` | Condor submit file (queue 1, temporary) |

### Histogram filling (Stage 2)

| File | Role |
|------|------|
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsimOverlay.cxx` | RDF filler: per-centrality filters, 6 quality categories |
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsim.cxx` | Parent: fullsim-specific hist vars (reco effcy, det response) |
| `RDFBasedHistFilling/RDFBasedHistFillingPythia.cxx` | Grandparent: pair categorization (op, single_b) |
| `RDFBasedHistFilling/var1D_pythia_fullsim.json` | 1D variable definitions for fullsim |

### Plotting (Stage 3)

| File | Role |
|------|------|
| `plotting_codes/reco_effcy/PythiaFullsimRecoEffPlotter.cxx` | Base + PP + Overlay plotter classes |
| `plotting_codes/reco_effcy/plot_reco_effcy_pythia_fullsim_pp24.cxx` | Standalone pp caller script |

### Shared

| File | Role |
|------|------|
| `MuonObjectsParamsAndHelpers/FullSimSampleType.h` | Enum + I/O path helpers |
| `MuonObjectsParamsAndHelpers/PbPbBaseClass.h` | Centrality binning (CRTP): `{0,5,10,20,30,50,80}` |
| `MuonObjectsParamsAndHelpers/MuonPairPythia.h` | `MuonPairPythiaFullSimOverlayWTruth`, `MuonPythiaFullSimOverlayWTruth` |

## Class hierarchy

```
PythiaFullSimOverlayAnalysis
  inherits PythiaAlgCoreT<MuonPairPythiaFullSimOverlayWTruth, MuonPythiaFullSimOverlayWTruth, Self,
                           PythiaFullSimExtras, PythiaFullSimOverlayExtras, PythiaTruthExtras>
  inherits PythiaFullSimExtras       -- reco-truth matching, medium/tight quality cuts
  inherits PythiaFullSimOverlayExtras -- overlay branches (centrality, FCal_Et, beam type)
  inherits PythiaTruthExtras          -- truth HF ancestry tracing

RDFBasedHistFillingPythiaFullsimOverlay
  inherits RDFBasedHistFillingPythiaFullsim  -- fullsim hist vars
  inherits PbPbBaseClass                      -- centrality bins

PythiaFullsimRecoEffPlotterOverlay
  inherits PythiaFullsimRecoEffPlotterBase   -- all plot logic
  inherits PbPbBaseClass                      -- centrality loop
```

## Pipeline

```bash
./pipelines/pipeline_pythia_fullsim_overlay.sh hijing [--dry-run]
```

Modes: `hijing`, `zmumu`, `data`. The `--dry-run` flag skips Condor submission and runs stages 2-5 using existing NTP output.

### Stage 1: NTuple processing

**What it does:** Reads NTUP files containing Pythia dimuon events overlaid on HIJING (or other) backgrounds. For each event, traces truth muon ancestry to tag HF origin, matches truth muons to reconstructed muons, applies quality selections, and writes flat `muon_pair_tree_sign1/sign2` trees.

**Overlay-specific processing:**
- Reads centrality (`avg_centrality`), calorimeter energy (`FCal_Et_A/C`), beam type
- Only processes pp beam type (`overlay_only_pp = true` by default)
- Weight: `ami_weight * nominal_beam_ratio / N_beam` (accounts for beam-luminosity ratio in overlay)
- Sets `fill_kn_trees_fullsim = true` to process all 6 kn ranges in a single job

**Condor:** Single job (`queue 1` in `.sub` file). Worker script: `run_pythia_fullsim_overlay_condor.sh <sample_type>`.

**Output:** `muon_pairs_pythia_fullsim_{label}_no_data_resonance_cuts.root`

### Stage 2: Histogram filling

**What it does:** Reads the flat trees via RDataFrame. Creates per-centrality filters on `avg_centrality` using bin edges `{0,5,10,20,30,50,80}`. For each centrality bin and pair category (`_op`, `_single_b`), fills histograms under 6 quality categories:

| Category suffix | Selection |
|----------------|-----------|
| (none) | All pairs in centrality bin |
| `_pass_medium` | `pair_pass_medium` |
| `_pass_tight` | `pair_pass_tight` |
| `_pass_signal_truth` | Truth-level signal region |
| `_pass_medium_and_signal_truth_and_reco` | Medium + truth+reco signal |
| `_pass_tight_and_signal_truth_and_reco` | Tight + truth+reco signal |

Signal region: `truth_minv in [1.08, 2.9]`, `truth_pair_pt > 8`, `truth_pair_eta < 2.2`, `truth_dr > 0.05`.

Also fills inclusive (no centrality cut) histograms via the parent `RDFBasedHistFillingPythiaFullsim`.

**Compilation:**
```bash
cd RDFBasedHistFilling
root -l -b -q '.L RDFBasedHistFillingPythiaFullsimOverlay.cxx+'
```

**Output:** `histograms_pythia_fullsim_{label}_no_data_resonance_cuts.root` (~2000 keys: 672 1D + 582 2D + 582 3D histograms + efficiency projection graphs)

### Stage 3: Plotting

**What it does:** Reads histograms, produces reco efficiency and detector response plots. Runs 4 configurations (medium/tight WP x with/without signal cuts), each covering inclusive + 6 centrality bins = 7 directories.

**Compilation:**
```bash
cd plotting_codes/reco_effcy
root -l -b -q '.L PythiaFullsimRecoEffPlotter.cxx+'
```

**Output structure:**
```
{data_dir}/plots/
  {plot_dir_prefix}_reco_effcy_plots/medium/   -- reco efficiency, all centrality bins
      reco_effcy_truth_pair_pt_ss_op_ctr0_5_mediumWP.png
      reco_effcy_truth_pair_pt_ss_op_ctr5_10_mediumWP.png
      ...
  {plot_dir_prefix}_reco_effcy_plots/tight/
  {plot_dir_prefix}_reco_effcy_plots_signal_cuts/medium/
  {plot_dir_prefix}_reco_effcy_plots_signal_cuts/tight/
  {plot_dir_prefix}_det_resp_plots/medium/     -- detector response, ctr_tag in filename
      pair_pt_truth_reco_compr_ctr0_5_mediumWP.png
  {plot_dir_prefix}_det_resp_plots/tight/
```
All centrality bins go into a single flat directory; the bin tag (`_ctrX_Y`) is embedded in each filename. No inclusive (centrality-integrated) plots are produced for the overlay.

### Validation between stages

The pipeline validates output at every transition:
1. After NTP: ROOT file integrity check (non-zombie, has keys) + `muon_pair_tree_sign1/sign2` non-empty
2. After RDF: ROOT file integrity + key count > 0
3. Plotting failures cause immediate exit

## Centrality binning

From `PbPbBaseClass.h`:

| Bin | Edge range | Tag string | Cross-section weight |
|-----|-----------|------------|---------------------|
| 0 | [0, 5) | `_ctr0_5` | `sigma_ctr[0]` |
| 1 | [5, 10) | `_ctr5_10` | `sigma_ctr[1]` |
| 2 | [10, 20) | `_ctr10_20` | `sigma_ctr[2]` |
| 3 | [20, 30) | `_ctr20_30` | `sigma_ctr[3]` |
| 4 | [30, 50) | `_ctr30_50` | `sigma_ctr[4]` |
| 5 | [50, 80) | `_ctr50_80` | `sigma_ctr[5]` |

Centrality is stored as `avg_centrality` (integer percentage) in the overlay NTUP. The HIJING test sample uses a 0–5 fm impact parameter selection (ultra-central), so only ctr [0,5) (~81%) and [5,10) (~18%) are populated; all peripheral bins are empty for this sample.

## Overlay truth pair structure and HF-ancestry

The overlay NTUP `truth_mupair_*` pairs are **Pythia×Pythia** (all barcodes < 200,000). The skimmer pairs only stable truth muons; HIJING hard-scatter particles are not stable muons.

**Three bugs were fixed in Apr 2026 (see [overlay_combinatorial_background.md](overlay_combinatorial_background.md) for details):**

1. **`truth_parents` CollectionProxy missing**: Fixed by `gInterpreter->GenerateDictionary("vector<vector<int>>","vector")` in `InitInputExtra`.

2. **Barcode collision in `GetParticleIndex` cache**: The overlay truth table has ~600–1200 duplicate barcodes per event (Pythia muon + Geant4 secondary assigned the same bc). The cache was built last-writer-wins, so lookups returned the wrong particle. Fixed with `emplace` (first-writer-wins).

3. **`m_eldest_bhadron_barcode` stored beam barcode**: In overlay the b-quark is absent; the chain ends at `B-meson → proton(bc=1)`. The B-meson barcode (not the proton's bc=1) is now stored for same-b matching.

**After fixes**, 500-event test (OS pairs, kn=0): `direct_b` 52%, `from_same_b` 30% — consistent with PP fullsim.

The residual ~15% `s_light` is genuine Geant4 combinatorial background (pion/kaon decay muons passing the pT > 4 GeV cut). These pairs appear in `_op` efficiency but not `_single_b`.

## Running manually (without pipeline)

```bash
source ~/setup.sh

# Stage 1: NTuple processing
cd NTupleProcessingCode
root -b -l <<'EOF'
.L PythiaAnalysisClasses.h
PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
py.fill_kn_trees_fullsim = true;
py.Run();
.q
EOF

# Stage 2: Histogram filling
cd ../RDFBasedHistFilling
root -l -b <<'EOF'
.L RDFBasedHistFillingPythiaFullsimOverlay.cxx+
{ RDFBasedHistFillingPythiaFullsimOverlay fs;
  fs.fullsim_sample_type = FullSimSampleType::hijing;
  fs.Run(); }
.q
EOF

# Stage 3: Plotting
cd ../plotting_codes/reco_effcy
root -l -b <<'EOF'
.L PythiaFullsimRecoEffPlotter.cxx+
{ gROOT->SetBatch(kTRUE);
  PythiaFullsimRecoEffPlotterOverlay pl(FullSimSampleType::hijing, false, false);
  pl.Run(); }
.q
EOF
```
