# Plotting Conventions

## Style

- ATLAS style is available at `Analysis/AtlasStyle.C` — load with `SetAtlasStyle()` from `AtlasStyle.h`
- Shared plot helpers: `Analysis/Utilities/PlotUtils.h` (SetStyle, LogAx, adjustLogXRange)
- Variable-to-title and log-scale maps: `Analysis/Utilities/PlotCommonConfig.h`

## Output format

- Save plots as PNG only (`c->SaveAs("....png")`)
- Do NOT also save PDF unless explicitly requested

## Cross-section plots

- PbPb cross-section plots are always combined (all years together), never per individual year

## Differential scaling

- For differential plots, use `h->Scale(NORM_FACTOR, "width")` — never manual bin-width division

## Labels & legends (physical wording — never arbitrary)

Every text drawn on a canvas (TLegend entries, axis/histogram titles, TLatex/TPaveText)
must use PHYSICAL wording, never internal/arbitrary code conventions:
- **Never `sign1`/`sign2`, never bare `OS`/`SS` as main text.** Write out "same sign" /
  "opposite sign" for dimuon pairs, or `#mu^{+}` / `#mu^{-}` for single muons. `OS`/`SS`
  are allowed ONLY as subscripts on a defined symbol (e.g. `G_{OS}`, `N_{SS}`,
  `k = G_{SS}/G_{OS}`).
- **Ranges use physical values with units, never bin numbers.** Label an interval as e.g.
  `8 < p_{T}^{pair} < 15 GeV`, never `bin 3` / `bins 1-5` / `slice 2`. Bin indices are
  arbitrary and change with the binning.
- **No raw tree/branch names or internal flags** (`df_op`, `_ss`, file suffixes) in labels —
  use the physics meaning.
- **Analysis-specific symbols carry their defining equation** on the plot (e.g.
  `k(m) = G_{SS}/G_{OS}`), unless the definition is very long. Domain-common quantities
  (`R_{AA}`, `v_n`, cross-section) need no on-plot definition.
- The nominal COARSE and fine-log pair-p_{T} binnings are read from `ParamsSet.h`
  (`pair_pt_coarse_bins` / `N_COARSE_PAIR_PT_BINS`), never re-invented per plot; label the
  resulting bins by their physical p_{T} ranges.

## MC per-event weights (MANDATORY)

Every MC-derived histogram (truth / fullsim / HIJING overlay / POWHEG) MUST be filled with the
per-event MC weight (generator `weight` / AMI xs·genFiltEff·beam_ratio/N / pT-hat-slice / dσ
weight), **regardless of normalization — including unit-area / shape-only / "normalized to unity"
overlays**. Never normalize a raw-count histogram: unweighted MC mixes pT-hat slices by generated
statistics, not cross-section, giving a distorted shape. Clone the WEIGHTED histogram, then scale
to unit area. (Enforced as criterion C6 in `physics-results-review.md`; memory
`feedback_hf_muons_not_prompt`.)

## Stacked histograms (THStack) — LINEAR y + magnitude ordering (MANDATORY)

A `THStack` shows composition as "parts of a whole": each component is the vertical THICKNESS
of its band and the reader adds bands to the total. Two hard rules preserve that reading:

- **LINEAR y-axis — never log.** On a linear axis equal vertical distance = equal amount, so a
  band's thickness is faithfully proportional to its contribution and its fraction of the total
  is readable by eye. On a **log** axis equal distance = equal *ratio*: a contribution Δ stacked
  on a running base B spans `log(1+Δ/B)`, which shrinks as B grows — so the same Δ looks tall at
  the bottom and invisible on top of a large base, the "sum of bands = total" reading breaks, and
  reordering the same components makes identical data look completely different. (Log y is for
  line/marker OVERLAYS comparing shapes — NOT for stacks.)
- **Order small/flat at the BOTTOM, most-abundant (largest integral) / steepest-peak at the TOP.**
  A band's bottom edge is the cumulative sum beneath it, so its drawn shape rides on that baseline.
  A small/flat component placed ON TOP of a tall steep peak is dragged into the peak's shape (looks
  peaked when it isn't). Put flat/small at the bottom (baseline ≈0 → true shape) and the big/steep
  peak on top (rides a nearly-flat base, undistorted; nothing above it to distort). Implement by
  sorting the `Add()` order by histogram **integral ascending** (smallest first = bottom); build
  the legend to match the visual top→bottom.

(Memory `feedback_thstack_linear_and_ordering`. This is the stack exception to
`feedback_log_scale_plots`, which otherwise wants log y for wide-dynamic-range 1D distributions.
Enforced as criterion C7 in `physics-results-review.md`.)
