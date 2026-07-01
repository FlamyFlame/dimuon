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
