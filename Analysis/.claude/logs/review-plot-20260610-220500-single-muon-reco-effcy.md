# Plot Review Log
**Task**: Review the single-muon reconstruction efficiency plots for both PP and HIJING overlay. PP plot: /usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/plots/pp24_single_muon_reco_effcy/single_muon_reco_effcy_vs_pt.png. Overlay plots (7 total): /usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/plots/hijing_overlay_pp24_single_muon_reco_effcy/single_muon_reco_effcy_vs_pt_ctr_inclusive.png, single_muon_reco_effcy_vs_pt_ctr0_5.png, single_muon_reco_effcy_vs_pt_ctr5_10.png, single_muon_reco_effcy_vs_pt_ctr10_20.png, single_muon_reco_effcy_vs_pt_ctr20_30.png, single_muon_reco_effcy_vs_pt_ctr30_50.png, single_muon_reco_effcy_vs_pt_ctr50_80.png. Source macro: plotting_codes/reco_effcy/plot_single_muon_reco_effcy.cxx. Physics: single-muon reco efficiency = pass_medium / all fiducial truth muons (pT > 4 GeV, |eta| < 2.4) vs truth pT in 9 coarse q*eta bins. For overlay: ev_centrality filter per centrality bin. PP: 431,424 entries, overlay: 107,796 entries (pp beam only).
**Log file**: review-plot-20260610-220500-single-muon-reco-effcy.md
**Started**: 2026-06-10T22:05:00Z
**Status**: APPROVED at iteration 1
**Iterations completed**: 1
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**:
None found.

**Numerical verification**:
| Quantity | Executor's value | Verified value | Source | Result |
|---|---|---|---|---|
| PP tree entries | 431,424 | 431,424 | muon_tree | MATCH |
| Overlay tree entries | 107,796 | 107,796 | muon_tree | MATCH |
| Centrality 0-5% | 87,855 | 87,855 | ev_centrality filter | MATCH |
| Centrality 5-10% | 19,218 | 19,218 | ev_centrality filter | MATCH |
| Centrality 10-20% | 723 | 723 | ev_centrality filter | MATCH |
| Centrality 20-30% | 0 | 0 | ev_centrality filter | MATCH |
| Centrality 30-50% | 0 | 0 | ev_centrality filter | MATCH |
| Centrality 50-80% | 0 | 0 | ev_centrality filter | MATCH |

**Summary**: All 8 plots (1 PP + 7 overlay) reviewed and approved. Efficiency values, axis labels, legends, binning, output directories, and physics all verified correct. Peripheral overlay bins (20-80%) correctly empty due to test sample centrality distribution.
