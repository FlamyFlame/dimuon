# Investigation Review Log
**Task**: pp24 data line in plots/mc_data_compr/signal/pair_pt_in_eta_subplots_mc_data_compr.png shows zero statistics in the 2nd-to-last pair-pT bin for pair-eta [-1.5,-1], but plots/single_b_analysis/pp24/pp_counts_pair_pt_in_eta_subplots.png (raw counts, same crossx binning) shows nonzero data count in that exact cell. Root cause needed; reco-eff correction cannot explain it.
**Log file**: review-investigation-20260903-132126-pp24-mc-data-zero-bin.md
**Started**: 2026-09-03T17:21:30Z
**Status**: IN PROGRESS
**Iterations completed**: 0
**Max iterations**: 5

## Iteration 1 (executor findings, pre-review)

**Hypothesis tested and CONFIRMED**: the "vanishing" bin is a pair-pT AXIS-BINNING MISMATCH
between the pp_counts plot family and the mc_data_compr plot family, not a data/selection/
efficiency bug. No pairs are actually lost anywhere.

**Files examined**:
- plotting_codes/mc_data_compr/plot_mc_data_pair_pt_in_eta.cxx (data hist = `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts`)
- plotting_codes/single_b_analysis/plot_pp_counts_pair_pt_in_eta.cxx (hist = `h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts`)
- plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx (confirms the DEFAULT/authoritative
  `pp24_crossx_pair_pt_in_eta_subplots.png` -- the one pp_counts' own docstring claims identical
  binning to -- is built from `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts`, i.e. pT_bins_120;
  the pT_bins_150 variant is opt-in (`use_pt_bins_150`, default false) and writes to a separate
  `_pt_150` output dir)
- RDFBasedHistFilling/RDFBasedHistFillingPP.cxx:754-897 -- both `h2d_counts_...` (line 835-837) and
  the DEFAULT `h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts` (line 809-811) are booked with
  `npt, ptbins` = **ParamsSet::pT_bins_120** (15 log bins, 8-120 GeV), on the SAME RDF node
  `df_single_b_crossx_weighted`. The MC-data-comparison-only histogram
  `h2d_crossx_pt_150_pair_eta_binned_w_signal_cuts` (line 895-897) is booked in a SEPARATE block
  with `npt150, ptbins150` = **ParamsSet::pT_bins_150** (15 log bins, 8-150 GeV) -- same bin COUNT
  (15) but a DIFFERENT edge array, extended so the axis can show the MC/POWHEG tail out to 150 GeV.
- MuonObjectsParamsAndHelpers/ParamsSet.h:629-630 -- `fillLogBinningArray(pT_bins_120, 15, 8, 120)`
  vs `fillLogBinningArray(pT_bins_150, 15, 8, 150)`; log10-based, confirmed via
  `fillLogBinningArray` body (lines 375-384).

**Quantitative check** (ROOT read of `dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root`,
pair-eta panel [-1.5,-1], the panel the user flagged):

pT_bins_120 (pp_counts / default pp24 crossx plot), raw counts:
  bin 13 [69.82,83.63]: N=20
  bin 14 [83.63,100.18]: N=10   <- the user's "nonzero second-to-last bin"
  bin 15 [100.18,120.00]: N=0

pT_bins_150 (`h2d_crossx_pt_150_...`, used ONLY by the mc_data_compr plot family), sum of weight:
  bin 12 [68.65,83.46]:  sumw=0.144
  bin 13 [83.46,101.47]: sumw=0.0843  <- physically the SAME range as pT_bins_120 bin 14, nonzero
  bin 14 [101.47,123.37]: sumw=0     <- the mc_data_compr "second-to-last bin" (15 total) -- a
                                          DIFFERENT, higher pT range that pT_bins_120 does not even
                                          resolve (falls inside pT_bins_120's own empty bin 15)
  bin 15 [123.37,150.00]: sumw=0

Both axes agree once compared by PHYSICAL RANGE, not by array index: 83-101 GeV is nonzero on both
axes; >101 GeV is zero on both (pT_bins_120 bin 15 [100.2,120] = 0 too). Because both axes happen to
have 15 bins, "second-to-last bin" LOOKS like the same comparison but is not: pT_bins_150's extra
30 GeV of headroom (120->150) shifts every edge upward, so its 14th bin covers 101.5-123.4 GeV
while pT_bins_120's 14th bin covers 83.6-100.2 GeV -- roughly a full bin's worth of shift by the
top of the spectrum where the log binning is coarsest.

**Root cause**: NOT a code/selection/efficiency bug and no data is dropped. It is a violation, in
spirit, of the CLAUDE.md BLOCKING binning rule ("1D/2D/3D views of the SAME quantity MUST use the
SAME binning" / never compare-by-position across differently-defined axes): the pp_counts plot and
the DEFAULT pp24 crossx plot both use `pT_bins_120`; the mc_data_compr plot family instead uses the
opt-in `pT_bins_150` variant exclusively (no pT_bins_120 option exists there), each with 15 bins but
different edges. The user's cross-plot bin-position comparison ("second-to-last of 15") therefore
compared two different physical pT ranges and read a genuine physical zero at 101.5-123.4 GeV as if
it were the same cell as the nonzero 83.6-100.2 GeV cell in pp_counts.

**Alternative explanations tested and ruled out**:
- Reco-eff correction turning nonzero->zero: ruled out per user's own correct reasoning (a
  multiplicative correction cannot zero out a nonzero raw entry), and confirmed structurally: the
  "second-to-last" pT_bins_150 bin is zero at the RAW weighted-count level (`sumw` before any
  extra correction beyond `crossx_weight_trig_corr`, which is already baked into `h2d_crossx_pt_150_...`)
  -- there is no separate stage that could zero a nonzero cell.
- Sign convention: not implicated -- both histograms are booked from the same OS node
  `df_single_b_crossx_weighted`.
- Log-scale rendering hiding a nonzero bar: ruled out -- verified directly from bin CONTENT, not
  the rendered PNG; the pT_bins_150 bin is a true zero in the source histogram.
- Panel/eta mismatch: ruled out -- both histograms' eta[-1.5,-1] panel resolves to bins 10-14 on
  their respective (agreeing, both `pair_eta_crossx_bins`, 48-bin) eta axis.

**No fix applied**: this is a plot-provenance/labeling finding, not a code defect to patch --
pT_bins_150 is a deliberate, documented, differently-scoped axis (comment at
RDFBasedHistFillingPP.cxx:770-773 and plot_mc_data_pair_pt_in_eta.cxx:57-66 both state the axis
source by name). Flagging for the user rather than unilaterally changing a binning per the
BLOCKING binning rule (never change/reinvent a binning without the user).
Reviewer agent spawned (a9c9556c84f4589e6), awaiting verdict.

## Iteration 1 review (reviewer verdict FAIL)
**Reviewer verdict**: FAIL
**Issues found**: 1 (WARNING)
1. [INVESTIGATION QUALITY / MINOR RIGOR GAP] The "genuinely zero" claim for pT_bins_150 bin 14
   [101.47,123.37] rests on an inference from the overlapping-but-not-identical pT_bins_120 bin 15
   [100.18,120.00], not a direct raw-count measurement, because no raw-count histogram exists on
   the pT_bins_150 axis. Fix: book a raw-count twin on the pT_bins_150 axis and rerun, or query the
   tree directly for a direct N in that exact bin.
**Numerical verification**: all cited numbers/edges MATCH (independently reproduced to <0.1%,
mostly exact). Only the inference gap above was flagged.

## Amendment (iteration 2)
Added `h2d_counts_pt_150_pair_eta_binned_w_signal_cuts` to
RDFBasedHistFilling/RDFBasedHistFillingPP.cxx (pure additive booking, same node
`df_single_b_crossx_weighted`, no weight column, on the pT_bins_150 axis) mirroring the existing
`h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts` (pT_bins_120) pattern. Reran
`RDFBasedHistFilling/run_crossx_hist_filling_pp24.sh` (recompiles RDFBasedHistFillingPP.cxx+ via
ACLiC, refills the full pp24 nominal crossx output) to materialize it and get a direct raw count
for pT_bins_150 bin 14 [101.47,123.37], eta[-1.5,-1] — the exact cell the reviewer flagged.

**Direct measurement result** (new `h2d_counts_pt_150_pair_eta_binned_w_signal_cuts`, read from the
freshly regenerated `dimuon_data/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal.root`,
mtime 2026-09-03 13:31), pair-eta [-1.5,-1] (y-bins 10-14), RAW pair counts on the pT_bins_150 axis:
  bin 12 [68.65,83.46]:   N=21
  bin 13 [83.46,101.48]:  N=10   <- matches pT_bins_120 bin 14's N=10 (same physical range)
  bin 14 [101.48,123.37]: N=0    <- DIRECT measurement, not inferred: the reviewer-flagged cell
  bin 15 [123.37,150.00]: N=0
(all-eta total for bin 14 = 16, i.e. genuinely sparse everywhere at this pT, not just this panel)

This is a direct raw count, not an inference from a neighboring axis -- closes the reviewer's
WARNING. Code change: additive booking of `h2d_counts_pt_150_pair_eta_binned_w_signal_cuts` in
RDFBasedHistFilling/RDFBasedHistFillingPP.cxx (pT_bins_150 unweighted twin of the existing
pT_bins_120 `h2d_counts_...`), then reran `run_crossx_hist_filling_pp24.sh` (recompile + full pp24
data crossx refill) to materialize it. No other histogram touched; `pp24_crossx_...` /
`pp_counts_...` outputs (pT_bins_120) and `h2d_crossx_pt_150_...` (the mc_data_compr plot's actual
data source) are numerically unchanged by this addition.

## Iteration 2
**Reviewer verdict**: PASS
**Issues found**: 0
**Details**:
None found. Code change independently confirmed via `git diff` (pure 8-line additive insertion,
no pre-existing line touched). New histogram correctly unweighted, on the pre-existing
`pT_bins_150` axis (no BLOCKING binning-rule violation -- axis was not newly invented). Output
file independently read; all bin-content numbers reproduced exactly. N=0 for the flagged cell
(pT_bins_150 bin 14, [101.475,123.374] GeV, eta[-1.5,-1]) is now a direct measurement, closing the
prior WARNING. Original root-cause conclusion (pT_bins_120 vs pT_bins_150 axis mismatch) confirmed
unchanged and correct.
**Numerical verification**: all cited numbers MATCH (git diff scope, bin contents 21/10/0/0/16,
cross-axis totals consistent, clean recompile, no RDF-exception-swallowing evidence).

**Status**: APPROVED at iteration 2
**Summary**: Root cause = pair-pT axis mismatch (ParamsSet::pT_bins_120 vs pT_bins_150, both 15
bins but different physical edges) between the pp_counts/default-crossx plot family and the
mc_data_compr plot family; NOT a data, selection, or efficiency bug. Confirmed by direct raw-count
measurement (new additive histogram) after an initial neighboring-axis-inference gap was closed.
