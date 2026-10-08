# Per-run signal-region yield / luminosity (abnormal-run QA)

## Objective
For pp24 and Pb+Pb 2023–26 data, plot (# signal-region events) / (integrated luminosity) per run to expose
abnormal runs. One PNG per dataset (5 total). Evidence only — **no run-exclusion recommendation**
(memory `feedback_run_exclusion_needs_experts`).

## Procedure (as run, 2026-10-08)
- Input: nominal `muon_pairs_*_mindR_0_02.root`, tree `muon_pair_tree_sign2` (= opposite sign); Pb+Pb NTuples
  already carry the 5-cut event selection (applied at NTuple stage, `PassEventSelHook`).
- Cuts (verbatim from `ParamsSet`, same as the crossx fill): `pair_pass_tight` (Tight WP),
  `SignalMinvCutExpr` (1.08–2.9), `SignalPairPtCutExpr` (>9), `FiducialGapCutExpr` on both muons,
  `PairFiducialEtaCutExpr` (|η^pair|<2.2); Pb+Pb also centrality acceptance 0 ≤ avg_centrality < 80.
- Yield: RAW, unweighted, counted as distinct events (run, LB, BCID, entry); pairs differ by <0.3 %.
- Luminosity: `Prescale Corrected` column of the analysis-trigger table in `IntNotes/data/luminosity/`
  (pp24 `HLT_2mu4_L12MU3V`, pb⁻¹; Pb+Pb `HLT_mu4`, µb⁻¹ → nb⁻¹ ÷1000).
- Not drawn: runs with zero lumi (461633) and `PbPbBadRuns()` (2023: 461674, 462964). Data and lumi run lists
  otherwise match exactly for all five datasets (pp24 8 runs; PbPb23 58; 24 14; 25 42; 26 35).

## Code / outputs
- `plotting_codes/data_qa/per_run_signal_yield_table.cxx` (RDF per-run counts), `plot_signal_yield_per_lumi_per_run.cxx` (plots); commit `c296228`.
- `~/usatlasdata/dimuon_data/plots/data_qa/signal_yield_per_lumi_per_run/` — PNGs + `*_yield_per_lumi.csv` (N, L, N/L, err, pull vs dataset average).
- Rerun: `root -l -b -q` a macro that `#include`s the plot .cxx and calls `plot_signal_yield_per_lumi_per_run()`.

## Results (|pull| > 4 vs dataset average N/L)
| Dataset | Run | N/L vs avg | Pull |
|---|---|---|---|
| pp24 | 488239 (N=699, L=0.26 pb⁻¹) | ≈2.1× | +20σ |
| PbPb23 | 463315 | ≈0.43× | −23σ |
| PbPb24, PbPb25 | none | | |
| PbPb26 | 522384 | ≈0.48× | −24σ |
| PbPb26 | 522408 | ≈0.50× | −30σ |
| PbPb26 | 522200 | ≈0.70× | −12σ |
| PbPb26 | 522298, 522628, 522925, 523138 | ≈1.07–1.09× | +4–5σ |

522384/522408 are the runs already covered by `pbpb2026_runs_522384_522408_ntrk_lowband.md`.
Causes of the others are NOT investigated.

## Remaining work
Hand the plots/CSVs to the experts; investigate 488239, 463315, 522200 only if requested.
`/review-plot` explicitly waived by the user (2026-10-08).
