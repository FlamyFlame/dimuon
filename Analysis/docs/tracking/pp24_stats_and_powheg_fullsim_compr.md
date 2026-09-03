# pp24 single-b signal STATISTICS + POWHEG FullSim in the MC/data comparison

Mode: **IMPLEMENTATION**. Opened 2026-09-03.

## Objective

Two deliverables requested by the user (2026-09-03):

**W1 — pp24 raw statistics of the single-b signal.**
1. `plots/single_b_analysis/pp24/pp_counts_pair_pt_in_eta_subplots.png` — RAW pp24 pair
   COUNTS in the single-b signal region, pair pT on the x axis, the 9 canonical pair-eta
   ranges as subplots, on **exactly the binning of the crossx plot**
   (`pp24_crossx_pair_pt_in_eta_subplots.png`).
2. `plots/single_b_analysis/pp24/pp_counts_pair_pt_in_eta.csv` — the same numbers as a table,
   pair pT in COLUMNS, pair eta in ROWS.

**W2 — POWHEG FullSim added to the single-b crossx comparison plot.**
3. `plots/mc_data_compr/signal/pair_pt_in_eta_subplots_mc_data_compr.png` (and its pair-eta
   integrated twin `pair_pt_mc_data_compr.png`, same macro, same 2D histogram) gains a
   **POWHEG FullSim pp17-conditions** curve, selected as truth single-b with the SAME
   data-derived signal cuts already applied to Pythia. pp17 conditions are used because no
   POWHEG FullSim pp24-conditions sample exists — **stated in the TLegend**.
4. Draw/legend order becomes **POWHEG, Pythia, pp24 data**, with the DATA drawn LAST so it is
   on top and never hidden by an MC curve.

Parent doc: `mc_data_compr_signal_generic_split.md` (DONE 2026-08-25) built the signal family
and is the authority on its cuts, binning, normalization and colour convention. This doc
extends it; nothing there is revised.

## Physics Procedure

### 1. Motivation

**W1.** Every number on the crossx plots is `sum 1/(eps_trig*eps_reco)/L`, so a bin can look
populated while resting on a handful of pairs. The statistical reach of the measurement — how
finely it can be binned, where the error bars are Poisson-dominated, where a bin must be merged
— is a property of the RAW COUNT, which no existing plot shows. This is the count.

**W2.** Pythia (LO 2->2 + parton shower) and POWHEG (NLO matrix element + shower) predict
different b-bbar kinematics, in particular the pair-pT spectrum, which is exactly the
observable of this figure. Putting both beside the data turns the comparison from a
"data vs one generator" statement into a generator-dependence statement.

### 2. Top-level statement

**W1.** Per (pair-pT, pair-eta) cell of the crossx axes:
```
N(pT, eta) = # of OS pairs surviving the pp24 single-b SIGNAL REGION, UNWEIGHTED
```
i.e. the same RDF node that fills `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts`
(`df_single_b_crossx_weighted`: Tight WP + `signal_cuts`) but filled with NO weight column, so
the bin content is an integer pair count and the error bar is sqrt(N). It is deliberately NOT
`h2d_..._corr_raw`, which is the count times `1/L_int` and carries a weighted error.

**W2.** Per (pair pT, pair eta) cell, on `ParamsSet::pT_bins_150` x `ParamsSet::pair_eta_crossx_bins`:
```
DATA   : dsigma/dX = (1/L)*SUM_{OS, signal region} 1/(eps_trig^pair * eps_reco^pair)   [pb]
PYTHIA : dsigma/dX = SUM_{truth single-b OS in the SAME region} w_AMI * 1000           [nb->pb]
POWHEG : dsigma/dX = SUM_{truth single-b OS in the SAME region} weight_norm            [already pb]
```
`weight_norm = (sigma_evt * eps_filt)/N_gen` is a cross-section in **pb**, established at source
in `mc_data_compr_signal_generic_split.md` Step 9 (mean `EventWeights[0]` ~ 1.1e7; read as nb
that would be an 11 mb b-bbar cross-section, 55x the total sigma_bb at 5.02 TeV). POWHEG
therefore receives normalization factor **1**, and must NOT get Pythia's nb->pb x1000.

### 3. Step-by-step method

**(a) The signal region** is the pp24 crossx selection, read from code and never retyped:
OS + Tight WP + `1.08 < m < 2.9` GeV + pair pT > 8 GeV + BOTH muons outside every window of
`ParamsSet::single_mu_fiducial_gap_cuts`. No dR cut. (Established from code in
`mc_data_compr_signal_generic_split.md` Step 1b.)

**(b) The POWHEG FullSim mirror** applies (a) in TRUTH quantities plus `from_same_b`, exactly
as the Pythia fullsim mirror does (`RDFBasedHistFillingPythiaFullsim.cxx:266,288-291`). Two
things in the existing POWHEG fullsim code do NOT match and are therefore NOT reused:
  * `pass_signal_truth` (the `.Define` in `CreateBaseRDFsPowhegFullsimExtra`) still applies the
    LEGACY one-sided `q*eta < 2.2`, retired on 2026-08-17 in favour of the three gap windows;
  * `df_single_b_weighted` (same function, just above that loop) is
    `from_same_b && truth_dr < 1.0`, whereas Pythia's is `from_same_b` alone.
  (Line numbers are deliberately NOT quoted: this task's own 80-line header comment shifted them
  once already, and a stale `:185` in a BLOCKING registry is worse than no number.)
  A NEW filter is added beside them; both are left byte-unchanged because they feed the POWHEG
  reco-efficiency / detector-response outputs. (Additive pattern, per parent doc D8.)
  The extra `truth_dr < 1.0` is expected to be INERT inside the signal region
  (m < 2.9 GeV at pair pT > 8 GeV forces dR <~ 2m/pT = 0.725) — to be VERIFIED numerically, not
  assumed.

**(c) Binning.** Nothing is invented or retyped. W1 uses the axes of the crossx histogram it
mirrors (`ParamsSet::pT_bins_120` x `ParamsSet::pair_eta_crossx_bins`, 9 panels from
`CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap`). W2's POWHEG histogram is booked
from the SAME named binnings as the Pythia and data partners (`pT_bins_150`,
`pair_eta_crossx`), so a bin-edge divergence is impossible and the plotter's `AssertSameAxes`
covers it.

### 4. Negative constraints

1. Do **NOT** width-scale or efficiency-correct the W1 counts — they are raw pair counts. The
   crossx plot beside them is the corrected object; that is the entire distinction.
2. Do **NOT** modify `RDFBasedHistFillingPowhegFullsim`'s existing `pass_signal_truth`,
   `df_single_b_weighted`, or any reco-eff / detector-response booking. Additive only.
3. Do **NOT** give POWHEG the nb->pb factor 1000 (see §2).
4. Do **NOT** present POWHEG FullSim pp17 as a pp24 prediction: it is sqrt(s_NN) = 5.02 TeV
   against the data's 5.36 TeV, and pp17 detector conditions. The legend says so.
5. Do **NOT** let an MC curve be drawn after the data — the data must be on top (user
   instruction).

## Implementation Plan

1. **W1a** `RDFBasedHistFillingPP.cxx`: book the UNWEIGHTED signal-region counts TH2 on the
   crossx axes. Per §2/§3c.
2. **W1b** rerun the pp24 data RDF crossx stage (`run_crossx_hist_filling_pp24.sh`).
3. **W1c** new macro `plotting_codes/single_b_analysis/plot_pp_counts_pair_pt_in_eta.cxx`
   -> PNG + CSV. Per §3c.
4. **W2a** `RDFBasedHistFillingPowhegFullsim.cxx` + `var1D_powheg_fullsim.json`: new
   gap-cut truth single-b signal filter + the shared-axis comparison histograms. Per §3b/§3c.
5. **W2b** rerun the POWHEG fullsim pp17 RDF stage.
6. **W2c** `plot_mc_data_pair_pt_in_eta.cxx`: POWHEG curve, draw/legend order, data on top,
   POWHEG/data ratio marker. Per §2/§4.
7. Reviews: `/review-analysis-code` (steps 1,4) and `/review-plot` (steps 3,6). Docs + commit.

## Design Decisions

**D1. POWHEG FullSim gets a NEW filter, `_single_b_pass_signal_truth_gapcut`, not an edit of the
existing one.** `RDFBasedHistFillingPowhegFullsim.cxx`'s `pass_signal_truth` (in
`CreateBaseRDFsPowhegFullsimExtra`; `:259` after this task's edits) still applies the one-sided
`q*eta < 2.2` retired on 2026-08-17, and its `df_single_b_weighted` (`:245-248`) is
`from_same_b && truth_dr < 1.0` whereas the Pythia FullSim partner is
`from_same_b` alone. Both feed the POWHEG reco-efficiency and detector-response outputs, so both
are left BYTE-UNCHANGED and the comparison gets its own node. Same additive pattern as parent-doc
D8. Measured: the dropped `truth_dr < 1.0` is exactly INERT in this region (517 459 signal pairs
with AND without it), as the kinematic bound dR <~ 2m/pT = 0.725 predicts.

**D2. POWHEG FullSim is normalized PER SAMPLE, with its own weight column.** The class-wide
`weight_norm = weight / SumMetaNentriesBeforeFilter(ALL input files)` divides by
`N_bb + N_cc = 9 785 325` although a cross-section is `sigma_bb + sigma_cc =
SUM_bb w/N_bb + SUM_cc w/N_cc`. That is an N-weighted AVERAGE, i.e. each production comes out
~2x under-normalized. It cancels in every ratio the class was written for (reco-eff numerator and
denominator carry the same weight), which is why it has never mattered; it does NOT cancel here,
and it does not even cancel approximately, because the single-b signal lives ENTIRELY in the bb
production (`from_same_b` needs a b-flavoured hadron ancestor, `PowhegTruthExtras.c:1088-1092`)
while cc supplies ~half the denominator. Measured: 517 459 signal pairs from bb, **0** from cc.
So the comparison histograms use `weight_norm_per_sample = weight / N_gen(this file)`, built with
RDataFrame's `DefinePerSample`; `weight_norm` and everything filled from it are untouched. Effect:
1768.22 pb -> **3527.25 pb**, a factor 1.995 = (N_bb+N_cc)/N_bb, exactly as predicted.
Recorded in `docs/powheg.md` §Weighting for any future absolute-normalization consumer.

**D3. POWHEG = GREEN in the signal family.** Fourth colour after data=red, Pythia single-b=black,
Pythia all-OS=blue, so that inside `plots/mc_data_compr/signal/` each colour means ONE thing
across every PNG. Deliberately NOT blue, which already means "Pythia, all OS" on the neighbouring
canvas. Markers also differ (POWHEG 22 triangle, Pythia 21 square, data 20 circle), so the curves
stay separable without colour. `McDataComprColors::kSignalPowheg`.

**D4. Data drawn LAST.** The frame is painted with `d->Draw("AXIS")`, then POWHEG, then Pythia,
then the data again with `"E,same"` — user instruction: the data must never be hidden by an MC
curve where the three cross. Legend order matches the draw order.

**D5. The counts histogram is a new UNWEIGHTED TH2, not a rescaled `_corr_raw`.** `cw_raw` is
`weight * (1/L_int)`, so the existing stage histogram is the count divided by the luminosity AND
carries a weighted error. Counting needs an unweighted fill: bin content integer, error sqrt(N).

**D6. The counts PNG reuses `SingleBCrossxPlotterBase::DrawPairPtByEta` with
`differential = false`.** Same projection code path, same panels, same styling as the crossx
figure it sits beside, so the two cannot describe different cells. The CSV repeats the projection
and then ASSERTS that the nine panels sum to the 2D integral, which is the property the 48-bin
pair-eta axis exists to guarantee.

## Progress Log

- 2026-09-03 Step 0: doc created after triage of INDEX.md + the full parent doc
  `mc_data_compr_signal_generic_split.md`.

- 2026-09-03 Step 1 DONE (W1a). `RDFBasedHistFillingPP.cxx` (after the crossx 2D booking):
  `h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts`, filled from the SAME node
  (`df_single_b_crossx_weighted`) on the SAME axes (`pms.pT_bins_120` x
  `ParamsSet::N_PAIR_ETA_CROSSX_BINS/MIN/MAX`) with NO weight column. Per D5.

- 2026-09-03 Step 2 DONE (W1b). pp24 data RDF crossx rerun (`run_crossx_hist_filling_pp24.sh`).
  New key present; integral **705 396**, entries 705 404. The crossx histogram beside it is
  unchanged at 5784.00 pb (`pT_bins_120` axis), so the new booking perturbed nothing.

- 2026-09-03 Step 3 DONE (W1c). New macro
  `plotting_codes/single_b_analysis/plot_pp_counts_pair_pt_in_eta.cxx` (class
  `PPSignalCountsPlotter : SingleBCrossxPlotterBase`), per D6. PNG via
  `DrawPairPtByEta(..., differential = false)`; CSV written by `WriteCsv`, which throws if the
  nine panels fail to sum to the 2D integral. Added to `pipeline_pp_crossx.sh` Stage 6 so it can
  never go stale relative to the crossx figure it mirrors.

- 2026-09-03 Step 4 DONE (W2a). `RDFBasedHistFillingPowhegFullsim.cxx`: file-local anonymous
  namespace (`McVsDataSignalFilter` = `_single_b_pass_signal_truth_gapcut`, `McVsDataVar1Ds`,
  `McVsDataVar2Ds`, `kNgenSampleColumn`, `kMcVsDataWeightColumn`); the new node in
  `CreateBaseRDFsPowhegFullsimExtra` (non-mixed only) with `DefinePerSample` normalization per
  D2; the fill in `FillHistogramsFullSim`, deliberately NOT inside the file's try/catch idiom, so
  a missing key throws instead of silently dropping the POWHEG curve. `var1D_powheg_fullsim.json`
  append-only, 21 -> 23 (`truth_pair_eta_crossx` -> `pair_eta_crossx`,
  `truth_pair_pt_log_150` -> `pT_bins_150`); no existing entry touched.
  **Pre-existing selections `pass_signal_truth` (`:185`) and `df_single_b_weighted` (`:171-173`)
  are byte-unchanged** (D1).

- 2026-09-03 Step 5 DONE (W2b). POWHEG fullsim pp17 RDF rerun. 396 -> **453 keys**; the run reads
  BOTH productions (`N_bb = 4 905 394`, `N_cc = 4 879 931`, printed by the new code). The three
  new histograms integrate to 3527.25 pb. First run (before D2) gave 1768.22 pb; the ratio is
  1.99486 = (N_bb+N_cc)/N_bb to 5 digits.
  Independent RDF cross-check on the raw trees: single-b signal pairs = **517 459 from bb, 0 from
  cc**, identical with and without `truth_dr < 1.0`, and `sum(weight)/N_bb = 3527.2` pb.

  **The on-disk file was STALE and the rerun refreshed it** (same shape as the parent doc's
  POWHEG-truth finding). `histograms_powheg_fullsim_pp17.root` was dated 2026-03-06 and predates
  several changes to the producer, so of the 57 new keys only 3 are this task's: the other **54
  come from the `truth_dr_4_0` variable added to `reco_effcy_var1Ds`/`var2Ds`** (1 x 18 filters +
  2 x 18 = 54, matching exactly), and the pre-existing `_pass_signal_truth` histograms were
  recomputed under the current cut set (commit `7b74dfc`, 2026-06-22, removed `truth_dr > 0.05`).
  **Blast radius: none.** Grep confirms the only consumer of the NON-mixed file is the new
  comparison macro; `PowhegFullsimRecoEffPlotter` reads the `_mixed*` outputs, which this run does
  not touch.

- 2026-09-03 Step 6 DONE (W2c). `McDataComprColors.h` gains `kSignalPowheg = kGreen+2` (D3);
  `McDataComprConfig.h` gains `PowhegFullsimFile()`; `plot_mc_data_pair_pt_in_eta.cxx` draws the
  third curve on both views, asserts its axes against the data's, orders POWHEG/Pythia/data with
  the data last (D4), and adds a second (green) ratio marker to every ratio pad. Legend fixes
  found by LOOKING at the rendered PNGs: the one-line POWHEG entry was clipped mid-word at
  nine-panel scale (now two lines), and the lower-left legend position ran through the falling
  spectrum in the outermost pair-eta panel (now upper right, where the 5-decade fall leaves the
  quadrant empty by construction).
  `plots/mc_data_compr/signal/README.md` updated: three-sample table, the pp17-conditions caveat
  as standing caveat 3, and the per-sample normalization statement.

## Results & Observations

### W1 — pp24 signal-region statistics

**705 396 OS pairs** in the single-b signal region (Tight WP, `1.08 < m < 2.9` GeV,
pair pT > 8 GeV, both muons outside the gap windows), against a cross-section of 5784.00 pb on
the same `pT_bins_120` axis (5784.06 pb on `pT_bins_150`, which keeps the > 120 GeV tail).
Entries 705 404, the 8-pair difference being the > 120 GeV overflow.

Distribution (from `pp_counts_pair_pt_in_eta.csv`):
* pair pT: 65 564 / 160 637 / 168 490 in the first three bins, falling to **17 pairs in the last
  bin (100-120 GeV) summed over ALL nine pair-eta panels** — 60 in 83-100 GeV, 175 in 69-84 GeV.
  So above ~70 GeV the measurement is single- to double-digit counts per pT bin even before the
  pair-eta split.
* pair eta: 130 039 in the central `[-0.5, 0.5]` panel down to 31 351 in `[2.0, 2.4]`.
* the emptiest CELL of the 9 x 15 grid is 0 (`[-2.4,-2.0]` and `[-1.5,-1.0]` at 100-120 GeV);
  several outer-eta cells at 83-100 GeV hold 2-3 pairs.

### W2 — POWHEG FullSim pp17 in the signal comparison

| | integrated OS cross-section | / data |
|---|---|---|
| pp data 2024 | 5784.06 pb | 1 |
| Pythia FullSim, single-b (pp24 cond.) | 7843.57 pb | 1.356 |
| POWHEG FullSim, single-b (**pp17 cond., 5.02 TeV**) | **3527.25 pb** | **0.610** |

Shape: the POWHEG/data ratio RISES monotonically with pair pT, 0.51 in the first bin to ~0.85 by
50 GeV and consistent with 1 above ~60 GeV, i.e. POWHEG is softer than the data. Pythia/data is
1.27-1.46 over the whole range (1.310 and 1.265 in the first two pT bins, ~1.45 from 15 to
60 GeV), i.e. flat within errors above ~15 GeV. Both trends are reproduced in all nine pair-eta
panels.

**A pair-eta trend that is probably the DATA's, not a generator's** (raised by /review-plot).
BOTH MC/data ratios rise together in the two outermost |eta| panels: Pythia 1.53 / 1.51 at
|eta| > 2 against 1.31-1.40 elsewhere, POWHEG 0.634 / 0.614 against 0.56-0.66. A rise shared by
two independent generators is evidence about the denominator, not the numerators — the natural
suspect is the data side's eps_trig / eps_reco at |eta| > 2, which is where the gap cut bites
hardest and where the turn-on fits and the reco-eff map are most extrapolated. Not investigated
here; recorded so it is not read as a generator statement.
Reading these numbers requires the caveats in `plots/mc_data_compr/signal/README.md`: the data
is a background-unsubtracted, non-unfolded reco-level OS yield, both MC curves are pure truth
single-b, and POWHEG is a 5.02 TeV pp17-conditions sample against 5.36 TeV data.

Counts-histogram cross-checks that passed (verified directly on the file):
* its axes are bit-identical to `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts` (15 x 48);
* the set of non-empty cells is IDENTICAL to the crossx histogram's — 0 cells where one has
  content and the other does not, which is what proves the two were filled from the same node;
* every bin content is an exact integer and every error is exactly sqrt(N);
* counts / `_corr_raw` = 400.412 (= 1/L_int) in EVERY bin (verified by the code reviewer), i.e.
  the two views of the same cells differ by exactly the luminosity factor and nothing else.

**One way in which N and dsigma do NOT describe the same pair set** (raised by
/review-analysis-code, and deliberate): a pair whose trigger efficiency cannot be evaluated gets
the sentinel `w_trig = 0` (`RDFBasedHistFillingPP.cxx:392-393`; the ACTIVE
`pp_trig_eff_highpt_jump.md` owns that pathology), so it is COUNTED here but contributes exactly
0 to every crossx histogram. N is the raw statistical reach of the SELECTION; dsigma is built
from the subset that could be corrected. Measured: **no cell is emptied by this** — 0 cells with
counts > 0 and crossx == 0 — so the effect is partial, never gross. Stated in the code, in the
CSV provenance line and here.

Cross-checks that passed:
* all three 2D histograms are edge-identical on both axes (`AssertSameAxes`, relative 1e-9);
* the POWHEG 2D integral equals its own 1D projections (3527.2494 vs 3527.2494 / 3527.2643, the
  latter being the pair-eta axis keeping the > 150 GeV tail the pT axis overflows);
* per-sample normalization verified by hand: `sum(weight) = 1.73026e10` over the bb file alone,
  / `N_bb = 4 905 394` = 3527.2 pb.

### Review findings and their disposition (2026-09-03)

Both mandated reviews ran as read-only subagents. `/review-analysis-code` -> **PASS-WITH-COMMENTS**
(it re-derived N_gen, the 517 459 / 0 bb/cc split, `Sum(weight)/N_bb = 3527.26 pb` and the
counts/`_corr_raw` = 400.412 identity with its own RDF job, and judged the D2 normalization fix
correct); `/review-plot` -> **PASS-WITH-COMMENTS**, no CRITICAL, no regression on the parent doc's
P1-P11.

| # | sev | finding | disposition |
|---|---|---|---|
| C1 | MED | `DefinePerSample` returned the FIRST `id.Contains(path)` hit. `RSampleInfo::Contains` is a SUBSTRING test and the same directory holds `..._part1.root` / `..._backup_before_dCache.root`; one future rename making an input path a substring of another would mis-normalize the whole curve silently. | FIXED: all matches counted, `throw` unless exactly 1. |
| C2 | MED | The unrecognised-sample `throw` fires INSIDE the RDF event loop, which ROOT swallows before exiting 0 — a weaker guard than it reads. | FIXED: new `AssertMcVsDataHistsFilled()`, called from `WriteOutputExtra` AFTER the loop; it throws if the 2D key is absent or empty and otherwise prints entries + integral. Confirmed firing: `517459 entries, integral 3527.25 pb`. |
| C3 | MED | The "not wrapped in try/catch so it throws" comment over-claimed: a bad WEIGHT COLUMN is caught and `continue`d inside `FillHistogramsSingleDataFrame`. | FIXED: comment narrowed to what `map_at_checked`/`Var1DSearch` actually guarantee, and points at C2's post-loop check for the rest. |
| C4/C5 | LOW | Stale line references (`:185`, `:171-174`, `:129/:131`) in the new comments, this doc, and `signal_selection_change_impact.md`. | FIXED: replaced by FUNCTION names in the code and in the BLOCKING registry, since this task's own header comment shifted them once already. |
| C6 | LOW | `h2d_counts_*` counts `w_trig = 0` sentinel pairs that contribute 0 to the crossx figure beside it. | FIXED as documentation (it is correct behaviour): stated in `RDFBasedHistFillingPP.cxx`, in the CSV provenance line, and in Results above. Measured: 0 cells emptied. |
| C7 | LOW | "matches Pythia BIT FOR BIT" is true of the kinematic+gap selection but `from_same_b` is a DIFFERENT upstream algorithm in the two generators. | FIXED: comment qualified, and added as standing caveat 4 of `plots/mc_data_compr/signal/README.md` — it matters because the figure's subject is generator dependence. |
| C8 | LOW | Pre-existing stale comments calling the pair-eta axis "44 uniform bins" / "45 edges" after the 2026-08-25 44->48 change — under the BLOCKING binning rule, a wrong comment about a binning is the exact failure mode. | FIXED at source in `RDFBasedHistFillingBaseClass.cxx` and `ParamsSet.h`. |
| P1 | MINOR | Two green POWHEG ratio points (~4.1) fell outside the common ratio frame and rendered as bare error-bar stubs — indistinguishable from a missing point. | FIXED: new `McDataComprRatio::DrawOutOfRangeMarkers()` paints a hollow triangle at the frame edge pointing the way the point went; called on both ratio series of both canvases. |
| P2 | MINOR | The new counts plot set exposed NO Medium/Tight WP config var (repo rule), and hard-coded `"PP 2024"` although `run_year` is a parameter. | FIXED: `wp` argument (`"tight"` default) switching the DATA file to `_medium_wp` with an actionable error if it does not exist; info line built from `run_year` + WP. Registered in `docs/muon_wp_registry.md` §3b. |
| P3 | MINOR | `SingleBCrossxPlotterBase::DrawPairPtByEta`'s legend ended at x2 = 0.93 NDC while the frame ends at 0.90 (default right margin), so the last character of every info line sat OUTSIDE the frame with the frame line through it — in the counts figure AND in the crossx figure. | FIXED in the base class: `x2 = 1.0 - gPad->GetRightMargin()`. `plots/single_b_analysis/pp24{,_pt_150}/` regenerated. **Pb+Pb crossx plots keep the old placement until their next replot** — cosmetic only, no content change. |
| P4 | MINOR | The legend named POWHEG's 5.02 TeV but left the data's energy implicit, so the contrast the caveat exists for had to be known already. | FIXED: `pp data 2024 (5.36 TeV)` and `Pythia, single-b (pp24 cond.)`. |
| P5 | LOW | The five other `signal/` PNGs predated their own (rewritten) data input. Content verified unchanged. | FIXED: `plot_mc_data_compr_signal.cxx+()` rerun; all seven now carry one timestamp. |
| P6 | LOW | The CSV carried 9x15 bare integers with no statement of what was counted. | FIXED: a `#`-prefixed provenance line naming sample, WP, the full selection, the "unweighted / not corrected / not subtracted" status and both binning sources. |
| P7 | INFO | Both MC/data ratios rise together at \|eta\| > 2 — evidence about the data's efficiencies, not about either generator; and "Pythia/data flat at 1.35-1.45" understated the first two pT bins. | RECORDED in Results above, with the corrected 1.27-1.46 range. |
| P8 | INFO | 3 of 394 positive points fall below the 9-panel common y frame (percentile range; printed to stdout). | ACCEPTED, unchanged from the parent design and auditable. |

## Remaining Work

1. **`weight_norm` in `RDFBasedHistFillingPowheg` remains wrong for any multi-production run**
   (see D2). Only the new comparison histograms were fixed, additively. Every other POWHEG
   FullSim histogram from a bb+cc run is still normalized by the summed denominator and is
   therefore ~2x low in absolute terms. Harmless in the reco-efficiency ratios it is used for;
   fix it at source (or forbid absolute readings) before any other consumer quotes an absolute
   number. The POWHEG **truth** run is unaffected (bb only).
2. **`RDFBasedHistFillingPowhegFullsim`'s `pass_signal_truth` still carries the retired one-sided
   `q*eta < 2.2`** (`:185/:186`), as does `RDFBasedHistFillingPowhegTruth.cxx:224-226`. Left
   byte-unchanged on purpose (D1) — they feed reco-eff / NLO-template inputs — but they are now
   two different definitions of "the signal region" living in one file. Registered in
   `docs/signal_selection_change_impact.md`.
3. **POWHEG FullSim exists only in pp17 conditions.** 5.02 TeV and a Run-2 detector against
   5.36 TeV data. Producing a pp24-conditions POWHEG FullSim would remove the caveat; until then
   the ratio 0.610 mixes a generator difference with a beam-energy difference.
4. **The counts figure has no SS twin.** The signal-region SS pair count would give the
   combinatorial scale of the same cells; `book_signal_region_1d` already fills SS 1D
   cross-sections, so only the unweighted 2D is missing. Not requested.

## Latest Stage

**DONE 2026-09-03.** All seven plan steps complete; both mandated reviews ran and every finding is
FIXED or explicitly accepted (table above); everything affected was recompiled and regenerated
after the fixes (POWHEG fullsim RDF, the counts figure, the pp24 crossx plot set, and all seven
`signal/` PNGs). Deliverables on disk:
* `plots/single_b_analysis/pp24/pp_counts_pair_pt_in_eta_subplots.png` + `.csv`
* `plots/mc_data_compr/signal/pair_pt_mc_data_compr.png` and
  `pair_pt_in_eta_subplots_mc_data_compr.png`, both with the POWHEG curve
