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
  resulting bins by their physical p_{T} ranges. **1D, 2D and 3D views of the same quantity
  MUST share one binning** — never build a 1D "slices" panel by regrouping a different, finer
  axis (see `.claude/CLAUDE.md` §Binnings).

## Publication standard: NO illustrative text on a plot (MANDATORY — `/review-plot` P1)

**A plot is presented to a physics audience. Anything you would not want a physics professor
to read off the figure must not be drawn on it.** The audience cannot see the code, the
tracking doc, or any .txt output, so a plot must stand alone and must contain only what a
published figure would.

**FAIL the plot if the canvas carries any of:**
- **Explanatory / tutorial sentences** about method or intent — "All three: MC direct
  P(mu4 | reco mu), no tag-and-probe.", "Ratio pad: MC / data, both series.",
  "SAME reco tag family, differ ONLY by the HIJING overlay.", "Residual != 0 comes from fit
  quality, binning ...". The ratio pad is labelled by its own axis title; the series are
  identified by the legend. Rationale belongs in the tracking doc.
- **Pointers to files, code or internal artefacts** — "see fit_report.txt", "(see the macro)",
  a ROOT file or histogram name, a method/mode identifier that exists only in the source
  (`polyu_fixedRp`, `sf_closure`, `do_step4`). If a number is worth showing, **draw the
  number**, never a pointer to where it lives.
- **Drawing/implementation asides** in legend entries — "(points only)", "(+ fit, dashed)",
  "(read back)", "measured, / plateau".
- **Apologetic or hedging prose** — "no plateau in this cell -- no fit",
  "(chi2/ndf undefined -- it passes through every point)". State the fact tersely ("no fit")
  or draw nothing.

Legitimate on-canvas text: axis titles with units, the legend (physical series names), the
sample/working-point headline, physical bin ranges, defining equations, fitted parameter
values with uncertainties, chi2/ndf **as a number**, and a terse record of any point pushed
off scale by an axis cap (dropping data silently is worse).

## Fits: the EQUATION is mandatory (MANDATORY — `/review-plot` P2)

**Any fitted or interpolated curve drawn on a plot must be accompanied by its exact equation**,
placed with (and above) the fitted parameter values, so the audience knows what was fitted and
what each parameter means.

- **FAIL** on abstract or abbreviated identifiers standing in for the function: a legend entry
  or annotation saying `powerlaw_floatRp fit`, `polyu_fixedRp`, `erf+log`, `Fermi fit` with no
  formula. Internal method names must not appear on the canvas at all.
- The legend entry for the curve is simply **`fit`** (or `fit` plus a physical qualifier when
  several fits are overlaid, e.g. `fit, corrected MC`).
- Write the equation in ROOT LaTeX with the SAME parameter symbols used in the value list, and
  define any auxiliary variable, e.g.
  `f(#DeltaR) = 1 + u^{2}(a_{2} + a_{3}u + a_{4}u^{2})` with `u #equiv max(0, 1 - #DeltaR/R_{p})`.
- Bare parameter names (`A`, `n`, `p`, `R_{p}`, `a_{2}`) are meaningless without the equation —
  they are acceptable **only** next to it.
- Fixed parameters should be identifiable as fixed (no uncertainty, or marked), so a reader does
  not mistake a constraint for a measurement.

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

### The wide-dynamic-range EXCEPTION (log y is allowed, and sometimes required)

**LINEAR is the default and the requirement. But linear is not a suicide pact:** when the stacked
components span so many decades that a linear axis renders the plot *unreadable* — every band except
the largest collapses onto the axis and carries no information at all — a **log y-axis is
permitted**, because an unreadable plot communicates nothing, which is strictly worse than a
distorted-but-legible one.

Rule of thumb: **≳3 decades of dynamic range** across the stacked components or across the x-range.
The canonical case is a steeply falling spectrum, e.g. the Pythia fullsim dσ/dp_T-per-pT-hat-slice
crossx stacks (`plot_pythia_fullsim{,_overlay}_kn_pt_crossx.cxx`), which span ~5 decades: on a
linear axis only the lowest-p_T bins are visible and the entire high-p_T tail — the physics of
interest — is a flat line on zero.

If you take the exception you MUST:
- **Say so on the plot or in the code** — a one-line comment stating the dynamic range and that
  linear was rejected as unreadable. An undocumented `SetLogy()` on a stack is still a FAIL.
- **Keep the magnitude ordering rule** (below) — it is independent of the axis and still applies.
- **Understand what you gave up:** on log y the band thickness is no longer proportional to the
  contribution, so the "parts of a whole" reading is gone. If the POINT of the plot is the
  composition/fraction, use linear and restrict the x-range (or plot fractions/ratios instead) —
  do NOT reach for log to dodge a hard plot.

(Memory `feedback_thstack_linear_and_ordering`. Linear-stack is the exception to
`feedback_log_scale_plots`, which otherwise wants log y for wide-dynamic-range 1D distributions;
this section is the exception-to-the-exception. Enforced as criterion C7 in
`physics-results-review.md`.)

## Axis range, binning & ratio panels (enforced as R1–R3 in `/review-plot`)

**R1 — the axis range must fit the data.** The data must occupy most of the axis. Do not
draw an axis far beyond where the sample has entries (e.g. muon p_T to 100 GeV when the
pT-hat slice dies out at 20 GeV — the points end up squeezed into a corner), and do not
choose a range so narrow that a significant fraction of entries lands in the
underflow/overflow. Pick the range from where the entries actually are.

**R2 — log axis ⇔ log binning.** For a *variable* axis (x in 1D; x and y in 2D), if the
axis is drawn log (`SetLogx`/`SetLogy`) then that axis's **binning must be log-spaced**
(geometric edges, e.g. `edges[i] = lo * pow(hi/lo, i/nb)`), and vice versa. Uniform bins
under a log axis make bin widths visually deceptive. (This is the converse of
`feedback_log_scale_plots`, which requires a log axis for log-binned variables.)

**R3 — put the ratio on the plot when the ratio is the point.** If two/three curves share a
canvas and their ratio is informative — they are expected to agree within errors
(sample-vs-sample cross-checks, closure tests), or one is a reference for the others
(data/MC, corrected vs uncorrected, nominal vs variation) — the figure MUST carry a ratio
panel (bottom pad, dashed line at 1). Not required for curves merely overlaid for shape.
