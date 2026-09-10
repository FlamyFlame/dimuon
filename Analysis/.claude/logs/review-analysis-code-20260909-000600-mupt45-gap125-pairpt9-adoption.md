# Analysis Code Review Log
**Task**: muon pT 4.0->4.5 GeV + gap window (-1.30,-1.05)->(-1.25,-1.05) + signal pair-pT 8->9 GeV adoption; full Pb+Pb migration onto the fiducial + pair-level gap cuts and the 48-bin pair-eta axis; crossx default pair-pT axis swap to 16 log bins 9->150.
**Log file**: review-analysis-code-20260909-000600-mupt45-gap125-pairpt9-adoption.md
**Started**: 2026-09-09T00:06:00-04:00
**Status**: IN PROGRESS
**Iterations completed**: 2
**Max iterations**: 5

## Iteration 1
**Reviewer verdict**: FAIL
**Issues found**: 13 (1 CRITICAL, 9 WARNING, 3 INFO)
**Details**:
1. [BINNING / SILENT WRONG PHYSICS] CRITICAL — `RAA_plotting.cxx:45-49` groups the fine pair-pT axis
   with a hardcoded 15-bin map + retyped edges, but the histograms it reads are the UNSUFFIXED family
   re-booked onto `pT_bins_150` (16 bins, 9->150). Bin 16 (125.81-150 GeV) silently dropped from R_AA
   mode 3; the retyped edge list and all three drawn labels wrong. The two existing guards compare pp
   vs Pb+Pb and both axes moved together, so nothing fires.
2. WARNING — data fitter's three `pT_min` lead comments state "Low edge moved 4 -> 4.5", the opposite
   of D11 and of the code beneath them.
3. WARNING — data fitter Fermi form comment retypes the range as `[4.5,60]`; the data range is `[4,60]`.
4. WARNING — four live comments still name `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts`, which was
   renamed AND re-booked onto a different axis.
5. WARNING — closure coverage guard described as "Both are 8 -> 150 GeV"; both now start at 9.
6. WARNING — `MCTrigEffPairPtBinning.h` `Describe()` PRINTS "8-150 GeV" into job logs and the `_pt4bin`
   provenance.
7. WARNING — `PairTrigEffEvaluator.h:68-69` retypes `pair_pt > 8` (already fixed by the executor before
   the reviewer finished).
8. WARNING — two `RDFBasedHistFillingPP.cxx` comments still say the ntuple-inherent cut is `pt>4`.
9. WARNING — `RDFBasedHistFilling/README.md` describes `pT_bins_120` as "15 bins 8-120" and does not
   mention the new default.
10. WARNING — `docs/signal_selection_change_impact.md` §0 (designated ground truth) false in four ways.
11. INFO — duplicated near-contradictory comment in the Pb+Pb OS-SS block.
12. INFO — `pT_bins_80` is reused as a SINGLE-muon axis in `var1D_pythia_truth.json` (`pt1`); moving it
    to 9 makes its first edge the pair-pT cut. Pre-existing reuse.
13. INFO — trig-eff 1D `pair_pt_log` now starts at 9 while the trig-eff mode has no pair-pT cut, so
    pairs in [8,9) fall into underflow on that diagnostic. Display-only.

**Numerical verification**: ALL MATCH except issue 1. Reviewer independently compiled a read-only macro
against the edited ParamsSet.h and confirmed: pT_bins_150 = 16 bins 9->150; pT_bins_120 = 16 bins 9->120;
pair_pt_coarse_bins = 8 bins with edges 9, 12.7930, 18.1846, 25.8485, 36.7423, 52.2273, 74.2385,
105.5262, 150; nesting coarse[k]==fine150[2k] for all k; the 4-bin variant nests at 0,2,4,6,8;
pT_bins_8 = 20 bins 4.0->8.0 with 4.5 an EXACT interior edge at index 2 (D10 satisfied); pT_bins_40
18 bins 4.5->40; pT_bins_60 unchanged 8->60; pT_bins_80 12 bins 9->80;
single_mu_pt_coarse_bins {4.5,8,14,25,100}; pTbins[0]=4.5; signal_pair_pt_min=9; gap windows
(-1.25,-1.05)(-0.10,+0.06)(2.20,2.40); N_PAIR_ETA_CROSSX_BINS=48; pair_eta_fiducial_max=2.2.
Also confirmed `pair_pass_tight` is visible in GetColumnNames() so the new guard cannot false-throw.
MISMATCH: RAA_plotting fine-pT bin-count assumption 15 vs actual 16.

**Amendments applied (executor)**: all 1-11 fixed. Issue 1 fixed at root: `pair_pt_rebins` and
`pair_pt_labels` are now DERIVED from `ParamsSet::pT_bins_150` via new `MakePairPtGroups()` /
`MakePairPtLabels()` helpers (near-equal 3-way split, remainder to the low groups where the spectrum is
steepest), plus a throwing coverage assert in `HistProject()` that fires if the grouping does not tile
the histogram actually read from disk. Issue 6 fixed by COMPOSING the printed range from
`pms.pair_pt_coarse_bins.front()/.back()` so it can never drift again. Issues 12/13 left as INFO:
12 would be a binning change requiring user sign-off (BLOCKING Binnings rule) and is pre-existing;
13 is display-only. Also fixed, found while amending: the `RDFBasedHistFillingPowhegFullsim.cxx` file
header still claimed `pass_signal_truth` uses the legacy `q*eta < 2.2`, which D7 migrated.

## Iteration 2
**Reviewer verdict**: FAIL
**Issues found**: 12 (0 CRITICAL, 8 WARNING, 4 INFO)
**Details**: every WARNING is a stale claim in a comment, a doc, or a RUNTIME-EMITTED string — the
dominant defect class of this whole change. No code-correctness or physics defect was found.
1. W — `RAA_plotting.cxx:490` figure legend publishes `p_{T}^{pair}>8 GeV` on the final-results R_AA
   plot, 40 lines below code that reads `raa_pms.pair_pt_coarse_bins`. Binnings rule 5.
2. W — `DrCorrectionCrossxEvaluator.h:167` and `DrCorrectionCascadeEvaluator.h:220` PRINT
   "pair pT < 8 GeV" to stdout on every crossx run (+ matching comments at :138 / :178).
3. W — `RDFBasedHistFillingPythiaFullsim.cxx:23,46` still claim `pair pT > 8 GeV` and the derived
   bound 0.725 (correct value at 9 GeV is 0.644); the POWHEG sibling WAS fixed — a half-fix.
4. W — `FitMCSinglesEffcy.cxx:82-83` says it "mirrors the data fitter exactly" beside the one
   deliberate divergence, contradicting its own header and D11. Iteration-1 issue 2's defect class
   re-created in the MC file while fixing the data file.
5. W — `ParamsSet.h:72-74` still said Pb+Pb has NOT adopted `N_PAIR_ETA_CROSSX_BINS` and "still
   retypes 44, -2.4, 2.4" — false since D1, and inside the file that DEFINES the constant.
6. W — three more stale migration claims: `PairEtaPanelBins.h:28-33`, `RAA_plotting.cxx:417-436`
   (its throw message sends the reader after a code defect that no longer exists),
   `MCTrigEffPairSelection.h:137-138`.
7. W — D7 residue in POWHEG: `PowhegTruth.cxx:14,236-237` and `PowhegFullsim.cxx:45` still describe
   `pass_signal_truth` as the legacy `q*eta < 2.2`.
8. W — `docs/analysis_overview.md` §2 now contradicts `signal_selection_change_impact.md` §0,
   which this round rewrote: two designated ground-truth docs asserting opposite signal regions.
   Plus lower-impact stale claims in powheg.md, pythia_fullsim_overlay.md, muon_wp_registry.md and
   two `run_pythia_fullsim_*_mc_trig.sh` headers.
9-12. INFO — `ParamsSet.h:52` lead line still says "44 uniform bins"; `pT_bins_80` is a SINGLE-muon
   axis in two JSONs (broader than iteration-1 INFO 12); retyped reco-eff dR slices `{8,12,20,inf}`;
   and a run-configuration coupling (`trigger_effcy_calc` now carries a selection meaning via
   `pbpb_run3_mu4_force_nominal`) with no guard — all eight live PbPb scripts set it correctly.

**Numerical verification**: ALL MATCH. Reviewer recomputed every axis independently and traced the
new R_AA grouping by hand: base=2 rem=2 -> 3/3/2 coarse cells -> fine bins {1..6},{7..12},{13..16},
covered=16, maxbin=16, **top bin 125.81-150 GeV retained**; labels 9-25.8 / 25.8-74.2 / 74.2-150
match the real edges; the nesting assert provably fires on `pT_bins_120` (12.793 vs 12.44) and the
coverage assert on a 15-bin stale file; `ModePrepare()` precedes `HistProject()` and member-init
order is correct. Also verified: **16** `SignalPairPtCutExpr` call sites across 9 files with ZERO
retyped `pair_pt > 8` in live selections; 14 booking groups checked for `npt`/`ptbins` agreement
with **0 mismatches**; **0** live retyped 44-bin pair-eta axes in either spelling; **0** residual
muon-pT 4.0 in live selections.

**Amendments applied (executor)**: all 8 WARNINGs fixed, plus INFO 9. Emitted strings and legends
are now COMPOSED from ParamsSet rather than retyped. INFO 10/11/12 left deliberately: 10 and 11 are
pre-existing binning questions needing user sign-off under the BLOCKING Binnings rule; 12 is
verified-correct in all eight live scripts and is a guard-hardening suggestion, not a defect.
Three further stale `pair pT > 8 GeV` comments found by my own sweep in live plotting code
(`plot_mc_data_compr_signal.cxx`, `PlotMCDataComprBaseClass.h`, `pp24_secondary_vertex_stats.cxx`)
were fixed too. Tracking-doc hits were left alone on purpose: those are append-only historical
records of earlier rounds and are correct as history.

**Cross-session finding merged as R4**: the textual-vs-inherited rule and the resulting partition of
the six MC fill classes — trig-eff chain IMMUNE, the three fullsim classes (eps_reco, detector
response, template-fit MC, overlay per-centrality eps_reco) STALE until the NTuple rerun.
