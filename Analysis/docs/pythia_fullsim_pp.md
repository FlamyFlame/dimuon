# Pythia Fullsim PP-Condition Analysis

Pythia di-muon pairs from open heavy-flavour decays, fully simulated and reconstructed under pp (no overlay) conditions. Measures baseline muon-pair reconstruction efficiency and detector response without pile-up or heavy-ion background.

This serves as the reference for comparing overlay results (hijing, zmumu, data) to isolate the effect of the Pb+Pb environment on reconstruction.

## Key files

### NTuple processing (Stage 1)

| File | Role |
|------|------|
| `NTupleProcessingCode/PythiaAnalysisClasses.h` | Defines `PythiaFullSimAnalysis` |
| `NTupleProcessingCode/PythiaAlgCoreT.{h,c}` | Template algorithm core |
| `NTupleProcessingCode/PythiaFullSimExtras.{h,c}` | Reco-truth matching, medium/tight quality cuts |
| `NTupleProcessingCode/PythiaTruthExtras.{h,c}` | Truth ancestry tracing |
| `NTupleProcessingCode/run_pythia_fullsim.sh` | Standalone run script (not Condor) |

### Histogram filling (Stage 2)

| File | Role |
|------|------|
| `RDFBasedHistFilling/RDFBasedHistFillingPythiaFullsim.cxx` | RDF filler: reco effcy vars, detector response |
| `RDFBasedHistFilling/var1D_pythia_fullsim.json` | 1D variable definitions |

### Plotting (Stage 3)

| File | Role |
|------|------|
| `plotting_codes/reco_effcy/PythiaFullsimRecoEffPlotter.cxx` | `PythiaFullsimRecoEffPlotter` (pp child class) |
| `plotting_codes/reco_effcy/plot_reco_effcy_pythia_fullsim_pp24.cxx` | Standalone caller |

## Class hierarchy

```
PythiaFullSimAnalysis
  inherits PythiaAlgCoreT<MuonPairPythiaFullSimWTruth, MuonPythiaFullSimWTruth, Self,
                           PythiaFullSimExtras, PythiaTruthExtras>
  inherits PythiaFullSimExtras   -- reco-truth matching, quality cuts
  inherits PythiaTruthExtras      -- truth HF ancestry

RDFBasedHistFillingPythiaFullsim
  inherits RDFBasedHistFillingPythia  -- pair categorization

PythiaFullsimRecoEffPlotter
  inherits PythiaFullsimRecoEffPlotterBase  -- all plot logic
```

The pp fullsim has no centrality binning (no `PbPbBaseClass`). Histograms are inclusive only.

## Pipeline (manual, no automated script yet)

There is no automated pipeline script for pp fullsim. The test sample is small enough to run interactively. When full samples become available, a Condor-based pipeline analogous to the overlay one should be created.

### Stage 1: NTuple processing

```bash
source ~/setup.sh
cd NTupleProcessingCode
root -b -l <<'EOF'
.L PythiaAnalysisClasses.h
PythiaFullSimAnalysis py;
py.fill_kn_trees_fullsim = true;
py.Run();
.q
EOF
```

Or use the existing script:
```bash
./run_pythia_fullsim.sh
```

All 6 kn ranges are processed in a single job (no batching). The `fullsim_sample_type` defaults to `FullSimSampleType::pp`.

Output: `muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts.root` in `FullSimSampleInputDir(pp)`.

### Stage 2: Histogram filling

```bash
cd RDFBasedHistFilling
root -l -b <<'EOF'
.L RDFBasedHistFillingPythiaFullsim.cxx+
{ RDFBasedHistFillingPythiaFullsim fs;
  fs.Run(); }
.q
EOF
```

Output: `histograms_pythia_fullsim_pp24_no_data_resonance_cuts.root`

Fills reco efficiency histograms (truth vs reco, pass_medium/pass_tight) and detector response matrices for pair_pt, minv, dr.

### Stage 3: Plotting

```bash
cd plotting_codes/reco_effcy
root -l -b <<'EOF'
.L PythiaFullsimRecoEffPlotter.cxx+
{ gROOT->SetBatch(kTRUE);
  PythiaFullsimRecoEffPlotter pl(false, false);  // medium WP, no signal cuts
  pl.Run();
  PythiaFullsimRecoEffPlotter pl_tight(true, false);
  pl_tight.Run(); }
.q
EOF
```

Output: plots in `{data_dir}/plots/pp24_reco_effcy_plots/` and `pp24_det_resp_plots/`.

## Differences from overlay analysis

| Aspect | PP fullsim | Overlay fullsim |
|--------|-----------|-----------------|
| Background | None | HIJING / Zmumu / data |
| Centrality bins | No | Yes (6 bins) |
| Beam type filter | N/A | pp only (`overlay_only_pp`) |
| Weight | Standard MC weight | `ami_weight * nominal_beam_ratio / N_beam` |
| Extra branches | None | `avg_centrality`, `FCal_Et_A/C` |
| Plotter class | `PythiaFullsimRecoEffPlotter` | `PythiaFullsimRecoEffPlotterOverlay` |
