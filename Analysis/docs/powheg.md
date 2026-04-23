# Powheg Analysis (Truth and Fullsim)

## Overview

Powheg+Pythia8 NLO MC for bb and cc dimuon production. Used to derive reconstruction efficiency corrections and detector response matrices for the cross-section measurement.

Unlike Pythia (LO), Powheg is an NLO generator. This has several consequences throughout the analysis chain:

- **Separate bb/cc processing**: Powheg generates one heavy-quark pair per event. bb and cc samples are produced and processed independently via the `mc_mode` parameter, then combined at histogram filling or via the mixed-pair procedure.
- **Per-event cross-section weighting**: each event carries an NLO cross-section (`EventWeights[0]`) that must be multiplied by a generator-level filter efficiency and normalized by the total number of generated events.
- **Cross-section outlier cut**: NLO generators occasionally produce extreme-weight events that must be rejected.
- **Mixed-pair procedure**: because Powheg produces exactly one HQ pair per event, natural dimuon pairs are single-flavour only. To form bb+cc combined pairs, single muons are extracted and re-paired artificially.

## Sample variants

| Variant | Class | `mc_mode` | `run_year` | Purpose |
|---------|-------|-----------|------------|---------|
| Truth | `PowhegTruthAnalysis` | bb / cc | N/A | Generator-level kinematics |
| Fullsim (no truth) | `PowhegFullSimAnalysisNoTruth` | bb / cc | 17, 24 | Reco-only processing |
| Fullsim (with truth) | `PowhegFullSimAnalysisWTruth` | bb / cc | 17, 24 | Reco-truth matched |
| Fullsim overlay (no truth) | `PowhegFullSimOverlayAnalysisNoTruth` | bb / cc | 18, 23–26 | PbPb overlay |
| Fullsim overlay (with truth) | `PowhegFullSimOverlayAnalysisWTruth` | bb / cc | 18, 23–26 | PbPb overlay |

All constructors take `(file_batch, mc_mode, run_year, useLocal)`. Truth ignores `run_year`.

## Weighting

### Per-event weight formula

```
event_crossx = EventWeights[0]
event_weight = event_crossx * filter_effcy
```

where:
- `filter_effcy_bb = 0.003`
- `filter_effcy_cc = 0.001108`

### Normalization

At histogram filling (RDF stage), the per-event weight is normalized:

```
weight_norm = weight / nentries_before_cuts_sum
```

where `nentries_before_cuts_sum` is the sum of `nentries_before_cuts` from the `meta_tree` across all input ROOT files (bb + cc combined). This normalizes histogram integrals to cross-section × filter efficiency.

If `meta_tree` is unavailable (legacy files), falls back to unscaled `weight`.

### Cross-section outlier cut

Events with `|weight| > crossx_cut * filter_effcy` are rejected, where `crossx_cut = 5 × 10⁸`. This removes rare extreme-weight events produced by the NLO generator.

### Jacobian-corrected weight

For dR distributions, a Jacobian-corrected weight is defined:

```
weight_norm_over_truth_dr = (truth_dr > 0) ? weight_norm / truth_dr : 0.0
```

This compensates for the phase-space Jacobian when projecting from (deta, dphi) to dR near dR → 0. Applied to `dr` and `dr_zoomin` 1D histograms and to `(pair_pt, minv)` 2D histograms.

### Comparison with Pythia weighting

| Aspect | Powheg | Pythia truth | Pythia overlay |
|--------|--------|-------------|----------------|
| Per-event weight | `EventWeights[0] * filter_effcy` | `weight` | `ami_weight * nominal_beam_ratio / N_beam` |
| Normalization | `/ nentries_before_cuts_sum` | N/A | N/A |
| Outlier cut | Yes (`5 × 10⁸`) | No | No |
| Filter efficiency | Per-flavour (bb/cc) | N/A | N/A |

## Truth analysis

### Ancestor tracing

Powheg truth analysis performs detailed ancestry tracing for each muon pair, significantly more elaborate than Pythia's simple HF origin tagging.

**4 ancestor groups** (`AncestorGrouping`):
- `gg` — gluon-gluon fusion
- `qg` — quark-gluon scattering
- `single_g` — single gluon
- `qq` — quark-antiquark annihilation

**Additional truth classifications:**
- Oscillation tracking (B-meson mixing: 0, 1, or 2 oscillations)
- Same-parent analysis (whether both muons share a common ancestor)
- `from_same_b` flag: both muons from the same b-quark decay chain
- Parent grouping: direct b, b→c, direct c, s/light, single photon
- Hard-scattering analysis: 2→2 vs 2→3 processes

**QQPair output trees**: heavy-quark pair kinematic distributions (dR, dphi, minv, pT) categorized by sign, dphi direction, and ancestor group. Controlled by `output_QQpair_tree` flag.

### Truth RDF histogram categories

Three independent categorization schemes, each filled with both normal and Jacobian-corrected weights:

1. **General**: near/away × ss/op
2. **Origin-binned**: ss/op × (gg, gq, single_gluon, qqbar, incoming, single_b, all)
3. **Flavor-binned**: ss/op × (single_b, both_from_b, both_from_c, others)

Plus signal acceptance histograms: 2D (pT, eta) for numerator (signal cuts) and denominator (`from_same_b`).

### Truth pipeline (manual)

No automated pipeline script. Run manually per flavour:

```bash
source ~/setup.sh
cd NTupleProcessingCode

# bb truth (6 batches via Condor)
condor_submit run_powheg_truth_bb.sub   # queue 6

# cc truth (6 batches via Condor)
condor_submit run_powheg_truth_cc.sub   # queue 6
```

Each batch reads 5 input files (`file_batch * 5` files per batch). After all batches complete, hadd per-flavour files into combined outputs, then run RDF:

```bash
cd RDFBasedHistFilling
root -l -b <<'EOF'
.L RDFBasedHistFillingPowhegTruth.cxx+
{ RDFBasedHistFillingPowhegTruth runner;
  runner.Run(); }
.q
EOF
```

Output: `histograms_powheg_truth.root`

## Fullsim analysis

### NoTruth vs WTruth variants

Unlike Pythia fullsim (which always includes truth), Powheg fullsim has two variants:

| Variant | Class | `perform_truth` | Truth extras | Use case |
|---------|-------|-----------------|-------------|----------|
| NoTruth | `PowhegFullSimAnalysisNoTruth` | false | None | Reco-only efficiency |
| WTruth | `PowhegFullSimAnalysisWTruth` | true | `PowhegTruthExtras` | Reco-truth matched analysis |

WTruth inherits both `PowhegFullSimExtras` (reco-truth matching) and `PowhegTruthExtras` (ancestor tracing). NoTruth disables QQpair tree output and truth histograms via `PowhegFullSimExtras::InitParamsExtra()`.

### Single-muon pipeline

Extracts individual muon trees from fullsim output. Required as input for the mixed-pair procedure.

```bash
./pipelines/pipeline_powheg_fullsim_single_muon.sh
```

**Stages:**
1. Submit bb + cc Condor jobs (11 batches each): `run_powheg_fullsim_wtruth_bb_single_muon.sub`, `run_powheg_fullsim_wtruth_cc_single_muon.sub`
2. Wait for both clusters
3. Validate per-batch outputs (ROOT file integrity)
4. hadd per-batch files → combined single-muon trees per flavour
5. RDF histogram filling: `RDFBasedHistFillingPowhegFullsimSingleMuon`
6. Validate histogram output
7. Plot detector response (medium + tight WPs): `PowhegFullsimDetRespPlotterSingleMuon`

**Output files:**
- `single_muon_trees_powheg_{bb,cc}_fullsim_pp17.root` (combined per-flavour)
- `histograms_powheg_fullsim_single_muon_pp17.root` (histograms)

### Mixed-pair procedure

Because Powheg generates one HQ pair per event, natural muon pairs are single-flavour. To produce combined bb+cc pairs for efficiency studies, single muons from both flavours are extracted and re-paired artificially.

**Prerequisite**: single-muon pipeline must complete first.

```bash
./pipelines/pipeline_powheg_fullsim_mixed_pairs.sh [--smoke-test] [--mass-filter|--no-mass-filter] [--skip-mixing]
```

**Options:**
- `--smoke-test`: run 200-pair bb mixing test before full submission
- `--mass-filter` (default): restrict `truth_minv` to [1, 3] GeV; files in `mixed_mass_1_3GeV/`
- `--no-mass-filter`: no minv filter; files in `mixed/`
- `--skip-mixing`: skip Condor submission, validate existing files then proceed to RDF + plotting

**Stages:**
0. (Optional) Smoke test: 200-pair bb mixing, validate, clean up
1. Submit Condor mixing jobs: 480 batches via `run_powheg_fullsim_single_muon_mixing.sub`
2. Wait for cluster (with lone-job timeout handling)
3. Validate all 480 batch files + spot-check tree non-emptiness
4. RDF histogram filling: `RDFBasedHistFillingPowhegFullsim(17, true)` with `mixed_subdir`
5. Validate histogram output
6. Plot reco efficiency (medium + tight WPs): `PowhegFullsimRecoEffPlotter`

**Lone-job timeout**: if one Condor job remains running alone for > `LONE_JOB_TIMEOUT_SECONDS` (default 12 h) and its output file exists, the job is killed and the batch is added to `skip_batches_mixed`, which is passed to the RDF filler to skip that batch.

**`single_b` pair category**: in unmixed mode, opposite-sign pairs where both muons come from the same b-quark decay chain (`from_same_b && truth_dr < 1.0`). In mixed mode, `single_b` is not used as a pair category — only `_ss` and `_op` are available (because natural ancestry information is lost during mixing).

**Mixed vs unmixed differences in RDF filling:**

| Aspect | Unmixed | Mixed |
|--------|---------|-------|
| Pair categories | `_ss`, `_op`, `_single_b` | `_ss`, `_op` |
| Detector response | Yes | No |
| Input files | 2 (bb + cc hadded pair trees) | 480 batch files |
| Reco eff projections | `TGraphAsymmErrors` | `TH1D` (TH divide) |

### Fullsim pair processing (unmixed)

For non-mixed fullsim analysis (natural pairs), bb and cc are processed by separate Condor jobs then combined at the RDF stage:

```bash
# bb fullsim with truth (11 batches)
condor_submit run_powheg_fullsim_wtruth_bb.sub

# cc fullsim with truth (11 batches)
condor_submit run_powheg_fullsim_wtruth_cc.sub
```

After hadd, run RDF:
```bash
cd RDFBasedHistFilling
root -l -b <<'EOF'
.L RDFBasedHistFillingPowhegFullsim.cxx+
{ RDFBasedHistFillingPowhegFullsim fs;
  fs.Run(); }
.q
EOF
```

## Pipelines

| Pipeline | Script | Dependency | Stages |
|----------|--------|-----------|--------|
| Truth | (manual Condor + hadd + RDF) | None | NTP → hadd → RDF |
| Fullsim single-muon | `pipelines/pipeline_powheg_fullsim_single_muon.sh` | None | Condor NTP → validate → hadd → RDF → det resp plots |
| Fullsim mixed-pairs | `pipelines/pipeline_powheg_fullsim_mixed_pairs.sh` | Single-muon pipeline | (smoke test) → Condor mixing → validate → RDF → reco effcy plots |

**Pipeline dependency**: single-muon pipeline must complete before mixed-pairs pipeline can run. The combined single-muon trees are the input to the mixer.

## Key files

### NTuple processing

| File | Role |
|------|------|
| `NTupleProcessingCode/PowhegAnalysisClasses.h` | All 5 concrete Powheg classes |
| `NTupleProcessingCode/PowhegAlgCoreT.{h,c}` | Powheg algorithm core template |
| `NTupleProcessingCode/PowhegTruthExtras.{h,c}` | Truth ancestry tracing, ancestor grouping, QQPair |
| `NTupleProcessingCode/PowhegFullSimExtras.{h,c}` | Reco-truth matching, quality cuts |
| `NTupleProcessingCode/PowhegFullSimOverlayExtras.{h,c}` | Overlay-specific (centrality) |
| `NTupleProcessingCode/mix_powheg_single_muon_pairs.C` | Single-muon pair mixer |
| `NTupleProcessingCode/run_powheg_truth_{bb,cc}.{sh,sub}` | Truth Condor scripts |
| `NTupleProcessingCode/run_powheg_fullsim_wtruth_{bb,cc}.{sh,sub}` | Fullsim pair Condor scripts |
| `NTupleProcessingCode/run_powheg_fullsim_wtruth_{bb,cc}_single_muon.{sh,sub}` | Fullsim single-muon Condor scripts |
| `NTupleProcessingCode/run_powheg_fullsim_single_muon_mixing.{sh,sub}` | Mixing Condor scripts |

### Histogram filling

| File | Role |
|------|------|
| `RDFBasedHistFilling/RDFBasedHistFillingPowheg.{h,cxx}` | Base Powheg RDF class (I/O, weight_norm) |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegTruth.cxx` | Truth: origin/flavor/general histograms |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegFullsim.cxx` | Fullsim: reco effcy, det resp, mixed pairs |
| `RDFBasedHistFilling/RDFBasedHistFillingPowhegFullsimSingleMuon.cxx` | Single-muon det resp |
| `RDFBasedHistFilling/var1D_powheg_truth.json` | Truth 1D variable definitions |
| `RDFBasedHistFilling/var1D_powheg_fullsim.json` | Fullsim 1D variable definitions |
| `RDFBasedHistFilling/var1D_powheg_fullsim_single_muon.json` | Single-muon 1D variable definitions |

### Plotting

| File | Role |
|------|------|
| `plotting_codes/reco_effcy/PowhegFullsimRecoEffPlotter.cxx` | Pair reco efficiency plotter |
| `plotting_codes/reco_effcy/PowhegFullsimDetRespPlotterSingleMuon.cxx` | Single-muon det resp plotter |
| `plotting_codes/reco_effcy/plot_det_resp_N_reco_effcy_powheg_fullsim_pp17.cxx` | Standalone caller |
| `plotting_codes/reco_effcy/plot_det_resp_N_reco_effcy_powheg_fullsim_pp17_single_muon.cxx` | Single-muon standalone caller |

## Class hierarchy

```
PowhegTruthAnalysis
  inherits PowhegAlgCoreT<MuonPairPowhegTruth, MuonPowhegTruth, Self,
                           PowhegTruthExtras>
  inherits PowhegTruthExtras  -- truth ancestry, ancestor grouping, QQPair

PowhegFullSimAnalysisNoTruth
  inherits PowhegAlgCoreT<..., PowhegFullSimExtras>
  inherits PowhegFullSimExtras  -- reco-truth matching, quality cuts

PowhegFullSimAnalysisWTruth
  inherits PowhegAlgCoreT<..., PowhegFullSimExtras, PowhegTruthExtras>
  inherits PowhegFullSimExtras + PowhegTruthExtras

PowhegFullSimOverlayAnalysis{NoTruth,WTruth}
  same as above + PowhegFullSimOverlayExtras  -- centrality, overlay flags

RDFBasedHistFillingPowheg (base)
  ├── RDFBasedHistFillingPowhegTruth
  ├── RDFBasedHistFillingPowhegFullsim
  │     └── RDFBasedHistFillingPowhegFullsimOverlay
  └── RDFBasedHistFillingPowhegFullsimSingleMuon
```

## Differences from Pythia

| Aspect | Powheg | Pythia |
|--------|--------|--------|
| Generator order | NLO | LO |
| bb/cc processing | Separate (`mc_mode`) | Together |
| Event weight | `EventWeights[0] * filter_effcy` | `weight` or `ami_weight * beam_ratio` |
| Outlier cut | Yes (`crossx_cut = 5e8`) | No |
| Jacobian weight | `weight_norm_over_truth_dr` | `weight_norm_over_truth_dr` |
| Mixed pairs | Required (3-pipeline chain) | Not needed |
| Ancestor tracing | 4 groups + oscillation + QQPair | HF origin tag only |
| Fullsim variants | NoTruth + WTruth | WTruth only |
| Overlay variants | Yes (future) | Yes (hijing/zmumu/data) |
| File batching | 5 files per batch (truth: 6 batches, fullsim: 11 batches) | 1 kn range per batch (6 batches) |
