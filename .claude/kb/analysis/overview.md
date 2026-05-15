# Analysis Overview

## Physics goal

Low-mass dimuon continuum measurement in pp and Pb+Pb collisions at the LHC (ATLAS). Covers muon-pair production in the mass range ~1–3 GeV from open heavy-flavour (b/c-quark) decays.

## Signal and backgrounds

- **Signal**: opposite-sign muon pairs from correlated heavy-flavour decays (bb, cc)
- **Same-sign pairs**: combinatorial background; used for background estimation
- **Drell-Yan**: irreducible background at higher mass

## Datasets (Run 3, active)

| Label | Collision | Energy | Trigger mode | `trigger_mode` | File suffix |
|-------|-----------|--------|-------------|----------------|-------------|
| pp24 | pp reference | 5.36 TeV | 2mu4 | 3 | `_2mu4` |
| PbPb23 | Pb+Pb | 5.36 TeV | single mu4 | 1 | `_single_mu4` |
| PbPb24 | Pb+Pb | 5.36 TeV | single mu4 | 1 | `_single_mu4` |
| PbPb25 | Pb+Pb | 5.36 TeV | single mu4 | 1 | `_single_mu4` |

Run 2 datasets (pp17, PbPb15, PbPb18) exist for cross-checks but are no longer actively maintained.

## MC generators

| Generator | Type | Purpose |
|-----------|------|---------|
| Pythia8 (LO) | truth, fullsim pp, fullsim overlay | Baseline HF dimuon kinematics; reco efficiency via overlay |
| Powheg+Pythia8 (NLO) | truth, fullsim, fullsim overlay | NLO corrections; reco efficiency and detector response matrices |

## Three-stage analysis pattern

Every pipeline follows: (1) NTuple processing → flat trees, (2) RDF histogram filling, (3) Plotting.
See `Analysis/README.md` for full details.

## PbPb event selection

5 sequential cuts (each requires prior cuts to pass):
1. ZDC vs FCal banana cut (mu+5sigma per FCal slice, TGraph interpolation)
2. ZDC time box: |t_A| < cut_ns AND |t_C| < cut_ns
3. ZDC preamp sum: preamp_A < cut_A AND preamp_C < cut_C
   - PbPb23/24: hard-coded per-year symmetric scalars (300/300, 420/420 ADC)
   - PbPb25: per-run mu+7sigma Gaussian fit (42 runs, stored as TTree `t_preamp_per_run`)
4. nTrk HItight fraction vs total nTrk: lower bound (per-slice Gauss, mu-5sigma)
5. nTrk HItight vs FCal ET band (per-slice Gauss, mu+/-5sigma)

Derived cuts stored in `/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20YY/event_sel_cuts_pbpb_20YY.root`.

## PbPb centrality

- Run 3 (2023/24/25): all use PbPb2023 FCal-ET Glauber thresholds
- PbPb25: `centrality` branch is all zeros in skim; recalculated in NTuple processing via `UpdateCentrality()` using PbPb2023 thresholds
- PbPb24: uses `GetCentralityPbPb2024` in `UpdateCentrality()`
- Standard bins: {0, 5, 10, 20, 30, 50, 80}%

## FCal scaling

FCal ET reweighting procedure was fully removed. No per-event FCal weight exists.
`weight_for_RAA = CalculateWeightForRAA(avg_centrality, weight)` only.

## Cross-section plots

PbPb cross-section plots are always combined (all years together), never per individual year.
