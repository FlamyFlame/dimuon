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

## Pipeline: `pipelines/pipeline_pythia_fullsim_pp.sh` (2026-07-20)

End-to-end pp-fullsim pipeline for BOTH samples, selected by one argument:

```bash
# FULL production (the "_pdf" sample 803015-803020, read from the LOCALGROUPDISK farm):
./pipelines/pipeline_pythia_fullsim_pp.sh full
# TEST sample:
./pipelines/pipeline_pythia_fullsim_pp.sh test
# smoke test (tiny nevents_max, exercises every stage):
SMOKE_NEVENTS=2000 ./pipelines/pipeline_pythia_fullsim_pp.sh full --dry-run
# with the MC-based trigger efficiency (Stage 3 + Stage 10):
ENABLE_MC_TRIG_EFF=1 ./pipelines/pipeline_pythia_fullsim_pp.sh full
```

Stages: 0 preflight (FULL: all 6 farm slices + 6 AMI files, else FATAL) → 1–3 NTP (nominal,
single-muon, +MC-trigger variants if enabled) → 4 validate → 5 RDF hists → 6 validate → 7 reco-eff
+ det-response → 8 single-muon reco-eff + reco-distr → 9 crossx per pT-hat slice + kn table →
10 MC trig-eff (optional). Env: `USE_TIGHT_WP` (default 1=Tight), `ENABLE_MC_TRIG_EFF` (default 0),
`SKIP_NTP`/`SKIP_RDF`/`SKIP_PLOTS`.

**THE sample switch is `full`/`test`** — it drives the input dir, the AMI cross-section dir, the
isospin treatment (`FullSimSampleType.h`) AND the `_full` output suffix, so they cannot drift. The
FULL sample is pp-beam-only with isospin weight 1 (an HONEST pp cross-section); the TEST sample has
4 isospin beams (a production mistake) and its absolute σ carries the Pb 4:6:6:9 average — NOT a
physical pp σ, and the crossx plot says so. See `docs/tracking/pythia_fullsim_pp24_full_sample_skim.md`
and `docs/ami_weights.md`.

### Manual stages (the pipeline just orchestrates these)

### Stage 1: NTuple processing

```bash
source ~/setup.sh
cd NTupleProcessingCode
root -b -l <<'EOF'
.L PythiaAnalysisClasses.h
PythiaFullSimAnalysis py;
// ISOSPIN: pp CONDITIONS simulate pp collisions -> the pp beam alone, isospin weight 1.
// That is the DEFAULT. The pp24 TEST sample was produced with 4 isospin beams by MISTAKE,
// so a run over the TEST sample must opt back in explicitly:
py.setIsospinBeams(true);          // TEST sample only -- OMIT for the full sample
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
Within each `_reco_effcy_plots/{medium,tight}/` directory, plots are sorted into
category subdirs: `signed/` (1D SS/OS), `single_b_op_compr/`, `2d/`, `distr/`
(reco_distr macro), and `ranged/`. Det-resp plots stay in `_det_resp_plots/`
(not subdivided). See `docs/pythia_fullsim_overlay.md` for the full scheme.

## Differences from overlay analysis

| Aspect | PP fullsim | Overlay fullsim |
|--------|-----------|-----------------|
| Background | None | HIJING / Zmumu / data |
| Centrality bins | No | Yes (6 bins) |
| Beam type filter | N/A | pp only (`overlay_only_pp`) |
| Weight | Standard MC weight | `ami_weight * nominal_beam_ratio / N_beam` |
| Extra branches | None | `avg_centrality`, `FCal_Et_A/C` |
| Plotter class | `PythiaFullsimRecoEffPlotter` | `PythiaFullsimRecoEffPlotterOverlay` |
